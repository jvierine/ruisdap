#[cfg(test)]
mod product_tests {
    use super::*;
    #[test]
    fn reviewed_gain_is_recomputed_and_source_is_preserved() {
        let temp = tempfile::tempdir().unwrap();
        let web = temp.path().join("web");
        let cal = temp.path().join("calibration");
        let data = web.join("data/source/zenith");
        std::fs::create_dir_all(&data).unwrap();
        std::fs::create_dir_all(web.join("catalog")).unwrap();
        let t = 1712599200.;
        let mf = data.join("manifest.json");
        let manifest = json!({"date":"2024-04-08","antenna":"ZENITH","radar_frequency_hz":440.2e6,"range_m":[300000.,304500.],"fits":[{"id":"fit","transmitter_w":1e6,"start_unix_s":t,"end_unix_s":t+180.,"parameters":[[3500.,1000.,2000.,0.],[3500.,1000.,2000.,0.]]}],"parameter_names":["proxy","Ti","Te","v"],"raw_url":"raw.f32","plasma_url":"plasma.f32"});
        std::fs::write(&mf, serde_json::to_vec(&manifest).unwrap()).unwrap();
        let before = std::fs::read(&mf).unwrap();
        std::fs::write(web.join("millstone-config.json"),br#"{"datasets":[{"date":"2024-04-08","antenna":"ZENITH","data_root":"data/source/zenith/"},{"date":"2024-04-08","antenna":"MISA","data_root":"data/source/misa/"}]}"#).unwrap();
        let key = "a".repeat(64);
        std::fs::create_dir_all(cal.join(&key)).unwrap();
        let k = 4. * std::f64::consts::PI * 440.2e6 / 299792458.;
        let ne = 1e11;
        let frequency = ((ne * (1.602176634e-19_f64).powi(2)
            / (8.8541878188e-12 * 9.1093837139e-31)
            + 3. * k * k * 1.380649e-23 * 2000. / 9.1093837139e-31)
            .sqrt())
            / std::f64::consts::TAU;
        let notebook = json!({"revision":3,"points":[{"included":true,"antenna":"ZENITH","date":"2024-04-08","observation_id":"plasma","fit_id":"fit","gate":0,"range_m":300000.,"frequency_hz":frequency,"frequency_standard_error_hz":10000.,"multiplier":-123.,"plasma_density_m3":42.}]});
        std::fs::write(
            cal.join(&key).join("latest.json"),
            serde_json::to_vec(&notebook).unwrap(),
        )
        .unwrap();
        let mut index = Index::default();
        index.datasets.insert(
            "2024-04-08/ZENITH".into(),
            Dataset {
                raw: vec![],
                manifest: mf.clone(),
            },
        );
        index.plasma.insert(
            "plasma".into(),
            PlasmaSource {
                path: PathBuf::new(),
                key: String::new(),
                antenna: "ZENITH".into(),
                date: "2024-04-08".into(),
                start: t,
                end: t + 60.,
            },
        );
        let s = Store {
            index: Arc::new(index),
            cache: Arc::new(Mutex::new(VecDeque::new())),
            public: web.clone(),
            calibration_root: cal,
            token: Arc::new("test-token".into()),
        };
        let request = || CalibrationRequest {
            notebook_key: key.clone(),
            revision: 3,
            date: "2024-04-08".into(),
            antenna: "ZENITH".into(),
        };
        let expected = ion_proxy(ne, 2000., 1000., k) / 3500.;
        let preview = magic_constant_result(&s, &request()).unwrap();
        assert!((preview["magic_constant"].as_f64().unwrap() / expected - 1.).abs() < 1e-12);
        assert_eq!(preview["magic_constant_units"], MAGIC_CONSTANT_UNITS);
        assert_eq!(preview["independent_fit_count"], 1);
        assert_eq!(preview["point_count"], 1);
        assert_eq!(preview["points"][0]["transmitter_w"], 1e6);
        assert!(preview["between_fit_scatter"].is_null());
        assert!(!web.join("catalog/millstone-config.json").exists());
        assert_eq!(std::fs::read(&mf).unwrap(), before);
        assert_eq!(std::fs::read_dir(web.join("data")).unwrap().count(), 1);
        // Receiver signal and P_tx both double at fixed plasma. The imported
        // proxy has already divided out power, so C must remain unchanged.
        let mut twice_power = manifest.clone();
        twice_power["fits"][0]["transmitter_w"] = 2e6.into();
        let a0_k = 3500. * 1e6 / ((1. + 2000. / 1000.) * 300_f64.powi(2));
        twice_power["fits"][0]["parameters"][0][0] =
            (2. * a0_k * (1. + 2000. / 1000.) * 300_f64.powi(2) / 2e6).into();
        std::fs::write(&mf, serde_json::to_vec(&twice_power).unwrap()).unwrap();
        let doubled = magic_constant_result(&s, &request()).unwrap();
        assert!((doubled["magic_constant"].as_f64().unwrap() / expected - 1.).abs() < 1e-12);
        assert_eq!(doubled["points"][0]["transmitter_w"], 2e6);
        for invalid in [Value::Null, json!(0.), json!(-1.)] {
            let mut broken = manifest.clone();
            broken["fits"][0]["transmitter_w"] = invalid;
            std::fs::write(&mf, serde_json::to_vec(&broken).unwrap()).unwrap();
            assert!(magic_constant_result(&s, &request())
                .unwrap_err()
                .to_string()
                .contains("transmit power"));
            assert!(calibrate_density(&s, request()).is_err());
            assert!(!web.join("catalog/millstone-config.json").exists());
        }
        std::fs::write(&mf, &before).unwrap();
        let mut stale = request();
        stale.revision = 2;
        assert!(magic_constant_result(&s, &stale).is_err());
        let path = s.calibration_root.join(&key).join("latest.json");
        let mut duplicate = notebook.clone();
        duplicate["points"]
            .as_array_mut()
            .unwrap()
            .push(notebook["points"][0].clone());
        std::fs::write(&path, serde_json::to_vec(&duplicate).unwrap()).unwrap();
        let repeated = magic_constant_result(&s, &request()).unwrap();
        assert_eq!(repeated["magic_constant"], preview["magic_constant"]);
        assert_eq!(repeated["independent_fit_count"], 1);
        assert!(repeated["between_fit_scatter"].is_null());
        duplicate["points"][0]["included"] = false.into();
        duplicate["points"][1]["included"] = false.into();
        std::fs::write(&path, serde_json::to_vec(&duplicate).unwrap()).unwrap();
        assert!(magic_constant_result(&s, &request()).is_err());
        std::fs::write(&path, serde_json::to_vec(&notebook).unwrap()).unwrap();
        let result = calibrate_density(&s, request()).unwrap();
        assert_eq!(result["magic_constant"], preview["magic_constant"]);
        assert!((result["gain"].as_f64().unwrap() / expected - 1.).abs() < 1e-12);
        let config: Value = serde_json::from_slice(
            &std::fs::read(web.join("catalog/millstone-config.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(config["datasets"][1]["data_root"], "data/source/misa/");
        let path = web.join(config["datasets"][0]["data_root"].as_str().unwrap());
        let product: Value =
            serde_json::from_slice(&std::fs::read(path.join("manifest.json")).unwrap()).unwrap();
        assert!((product["fits"][0]["parameters"][0][0].as_f64().unwrap() / ne - 1.).abs() < 1e-12);
        assert_eq!(product["fits"][0]["parameters"][0][1], 1000.);
        assert_eq!(product["fits"][0]["uncalibrated_density_proxy"][0], 3500.);
        let h = hdf5::File::open(path.join("calibration.h5")).unwrap();
        assert_eq!(
            h.dataset("source_fit_parameters").unwrap().shape(),
            vec![1, 2, 4]
        );
        assert!(
            (h.dataset("calibrated_ne_m3")
                .unwrap()
                .read_raw::<f64>()
                .unwrap()[0]
                / ne
                - 1.)
                .abs()
                < 1e-12
        );
        assert_eq!(std::fs::read(mf).unwrap(), before);
        assert!(
            (h.dataset("magic_constant")
                .unwrap()
                .read_raw::<f64>()
                .unwrap()[0]
                / expected
                - 1.)
                .abs()
                < 1e-12
        );
        assert!(h.group("records/fit/fit").is_ok());
        assert_eq!(
            h.dataset("fit_transmitter_power_w")
                .unwrap()
                .read_raw::<f64>()
                .unwrap(),
            vec![1e6]
        );
        let meta = h
            .group("records/fit")
            .unwrap()
            .attr("metadata")
            .unwrap()
            .read_raw::<u8>()
            .unwrap();
        let meta: Value = serde_json::from_slice(&meta).unwrap();
        assert_eq!(meta["transmitter_w"], 1e6);
        assert_eq!(
            h.dataset("records/fit/fit/parameters").unwrap().shape(),
            vec![2, 4]
        );
        if let Ok(destination) = std::env::var("RUISDAP_TEST_EXPORT") {
            std::fs::copy(path.join("calibration.h5"), destination).unwrap();
        }

        let mut wrong = request();
        wrong.revision = 2;
        assert!(calibrate_density(&s, wrong).is_err());
    }
}
