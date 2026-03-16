# AutoQ-Rust: Automata-based quantum program verification
# Build: docker build "https://github.com/go2run/Autoq-Rust.git#claude/autoq-rust-analysis-R36fZ"
# Run:   docker run --rm -it <image>

FROM rust:1.85-bookworm AS builder

WORKDIR /app
COPY . .

RUN cargo build --release -p autoq-cli && \
    cargo test --release

FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends libgcc-s1 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/autoq /usr/local/bin/autoq

ENTRYPOINT ["autoq"]
