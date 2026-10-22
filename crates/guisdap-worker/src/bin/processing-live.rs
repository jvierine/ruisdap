#[path = "../live_monitor.rs"]
mod live_monitor;
#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let bind = std::env::var("RUISDAP_LIVE_BIND").unwrap_or_else(|_| "127.0.0.1:8093".into());
    let listener = tokio::net::TcpListener::bind(bind).await?;
    axum::serve(listener, live_monitor::router("/guisdap/api/live")).await?;
    Ok(())
}
