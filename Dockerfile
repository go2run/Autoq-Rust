# AutoQ-Rust Development Environment
# Build: docker build "https://github.com/go2run/Autoq-Rust.git#claude/autoq-rust-analysis-R36fZ" -t autoq-rust
# Run:   docker run --rm -it autoq-rust

FROM rust:1.85-bookworm

WORKDIR /home/user/Autoq-Rust

COPY . .

RUN cargo build --release -p autoq-cli && \
    cargo test --release

ENV PATH="/home/user/Autoq-Rust/target/release:${PATH}"

CMD ["bash"]
