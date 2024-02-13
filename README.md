# mat-vec-prod-exp
*Matrix vector product experiments*

Instructions:

- Just run `cargo test -- --nocapture` and the number of constraints will be printed
- Can comment & uncomment `lib.rs` lines `71` & `72` to test the different approaches:
	```rust
	let Az = mat_vec_mul_sparse_gadget(A, z);
	// let Az = handcrafted_A_by_z(cs, z)?;
	```
