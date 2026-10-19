"""A resumed analysis must not mix native numerical kernels or overwrite fits."""
import pathlib, subprocess, sys, tempfile, unittest
import h5py

REPO=pathlib.Path(__file__).resolve().parents[2]
class ResumeGuards(unittest.TestCase):
    def run_driver(self,root,*extra):
        return subprocess.run([sys.executable,str(REPO/'scripts/process-millstone-full-day.py'),
            '--root',str(root/'source'),'--date','2024-04-08','--output',str(root/'run'),
            '--historical',str(root/'historical'),'--publish',str(root/'public'),*extra],
            cwd=REPO,text=True,capture_output=True)
    def test_existing_products_require_explicit_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp);(root/'run/products').mkdir(parents=True)
            result=self.run_driver(root)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Use --resume-products explicitly',result.stderr)
    def test_changed_kernel_is_rejected_before_reading_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp);(root/'run/products').mkdir(parents=True)
            with h5py.File(root/'run/kernel.h5','w') as h:h.attrs['binary_sha256']='different-kernel'
            result=self.run_driver(root,'--resume-products')
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Native kernel differs',result.stderr)
            self.assertFalse((root/'run/manifest.json').exists())

if __name__=='__main__':unittest.main()
