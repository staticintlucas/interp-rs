# interp

[![Test status](https://img.shields.io/github/actions/workflow/status/staticintlucas/interp-rs/test.yml?branch=main&label=tests&style=flat-square)][tests]
[![Code coverage](https://img.shields.io/codecov/c/gh/staticintlucas/interp-rs?style=flat-square)][coverage]
[![Crate version](https://img.shields.io/crates/v/interp?style=flat-square)][version]
[![Rust version](https://img.shields.io/badge/rust-1.34%2B-informational?style=flat-square)][rust version]
[![Downloads](https://img.shields.io/crates/d/interp?style=flat-square)][downloads]

A Rust reimplementation of Matlab's `interp1` function for linear interpolation.

API documentation is available on [docs.rs][docs].

[tests]: https://github.com/staticintlucas/interp-rs/actions
[coverage]: https://codecov.io/gh/staticintlucas/interp-rs
[version]: https://crates.io/crates/interp
[rust version]: https://crates.io/crates/interp
[downloads]: https://crates.io/crates/interp
[docs]: https://docs.rs/stringslice/latest/stringslice/

## Usage

Add `interp` to your `Cargo.toml` file:

```toml
[dependencies]
interp = "0.1"
```

## Example

```rust
use interp::interp;

let x = vec![0.0, 0.2, 0.5, 0.8, 1.0];
let y = vec![0.0, 1.0, 3.0, 3.5, 4.0];

assert_eq!(interp(&x, &y, 0.35), 2.0);
```

Full API documentation is available on [docs.rs][docs].

[docs]: https://docs.rs/stringslice/latest/stringslice/

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

You can install the [pre-commit] hook (which checks formatting, etc) by running:

```bash
pip install -U pre-commit
pre-commit install
```

[pre-commit]: https://pre-commit.com/

## Licence

Licensed under either of

* Apache License, Version 2.0 ([LICENCE-APACHE] or [http://www.apache.org/licenses/LICENSE-2.0])
* MIT license ([LICENCE-MIT] or [http://opensource.org/licenses/MIT])

at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in
this crate by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without
any additional terms or conditions.

[LICENCE-APACHE]: LICENCE-APACHE
[http://www.apache.org/licenses/LICENSE-2.0]: http://www.apache.org/licenses/LICENSE-2.0
[LICENCE-MIT]: LICENCE-MIT
[http://opensource.org/licenses/MIT]: http://opensource.org/licenses/MIT
