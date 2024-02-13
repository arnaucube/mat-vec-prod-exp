#![allow(non_snake_case)]
#![allow(unused_doc_comments)]
#![allow(dead_code)]

use ark_ff::PrimeField;
use ark_r1cs_std::fields::nonnative::NonNativeFieldVar;
use ark_r1cs_std::{alloc::AllocVar, eq::EqGadget, fields::FieldVar, R1CSVar};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use core::marker::PhantomData;
use std::ops::Mul;

mod utils;
use utils::*;

/// - F stands for the field that we represent
/// - CF stands for the ConstraintField over which we do the operations

/// Implements the A * z matrix-vector-product by fixing the combinations of 'z'.
fn handcrafted_A_by_z<F: PrimeField, CF: PrimeField>(
    cs: ConstraintSystemRef<CF>,
    z: Vec<NonNativeFieldVar<F, CF>>,
) -> Result<Vec<NonNativeFieldVar<F, CF>>, SynthesisError> {
    let five = NonNativeFieldVar::<F, CF>::new_constant(cs.clone(), F::from(5u32))?;
    // directly hand-craft the output vector containing the operations in-place:
    Ok(vec![
        z[1].clone() + five.clone() * z[4].clone(),
        z[1].clone() + z[3].clone(),
        z[1].clone() + z[4].clone(),
        five * z[0].clone() + z[4].clone() + z[5].clone(),
    ]
    .clone())
}

/// Implements the A * z matrix-vector-product by doing the sparse matrix by vector algorithm, and
/// assuming that the elements of the matrix A are constants of the system.
pub fn mat_vec_mul_sparse_gadget<F: PrimeField, CF: PrimeField>(
    m: SparseMatrixVar<F, CF>,
    v: Vec<NonNativeFieldVar<F, CF>>,
) -> Vec<NonNativeFieldVar<F, CF>> {
    let mut res = vec![NonNativeFieldVar::<F, CF>::zero(); m.n_rows];
    for (row_i, row) in m.coeffs.iter().enumerate() {
        for (value, col_i) in row.iter() {
            if value.value().unwrap() == F::one() {
                res[row_i] += v[*col_i].clone(); // when value==1, no need to multiply by it
                continue;
            }
            res[row_i] += value.clone().mul(&v[*col_i].clone());
        }
    }
    res
}

/// Circuit that takes as constants the sparse matrix A, and as inputs the vectors z and y. It
/// computes the matrix by vector product between A and z, and checks that is equal to y
/// (ie. y == A*z)
struct MatrixVectorCircuit<F: PrimeField, CF: PrimeField> {
    _cf: PhantomData<CF>,
    pub A: SparseMatrix<F>,
    pub z: Vec<F>,
    pub y: Vec<F>,
}
impl<F: PrimeField, CF: PrimeField> ConstraintSynthesizer<CF> for MatrixVectorCircuit<F, CF> {
    fn generate_constraints(self, cs: ConstraintSystemRef<CF>) -> Result<(), SynthesisError> {
        // set A as circuit constants
        let A = SparseMatrixVar::<F, CF>::new_constant(cs.clone(), self.A)?;
        // set z and y as witness (private inputs)
        let z: Vec<NonNativeFieldVar<F, CF>> = Vec::new_witness(cs.clone(), || Ok(self.z.clone()))?;
        let y: Vec<NonNativeFieldVar<F, CF>> = Vec::new_witness(cs.clone(), || Ok(self.y.clone()))?;

        /// The next two lines are the ones that can be swapped to see the number of constraints
        /// taken by the two approaches:
        let Az = mat_vec_mul_sparse_gadget(A, z);
        // let Az = handcrafted_A_by_z(cs, z)?;

        Az.enforce_equal(&y)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_pallas::{Fq, Fr};
    use ark_relations::r1cs::ConstraintSystem;

    #[test]
    fn test_relaxed_r1cs_nonnative_matrix_vector_product() {
        let A = to_F_matrix::<Fq>(vec![
            vec![0, 1, 0, 0, 5, 0],
            vec![0, 1, 0, 1, 0, 0],
            vec![0, 1, 0, 0, 1, 0],
            vec![5, 0, 0, 0, 1, 1],
        ]);
        let z = to_F_vec(vec![1, 123, 35, 53, 80, 30]);
        let y = mat_vec_mul_sparse(&A, &z); // y = A*z
        println!("Matrix of size {} x {}", A.n_rows, A.n_cols);
        println!("Vector of size {}", z.len());

        println!(
            "Build the circuit that computes the matrix-vector-product over a non-native field"
        );
        let cs = ConstraintSystem::<Fr>::new_ref();
        let circuit = MatrixVectorCircuit::<Fq, Fr> {
            _cf: PhantomData,
            A,
            z,
            y,
        };
        circuit.generate_constraints(cs.clone()).unwrap();
        println!("Number of constraints: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap());
    }
}
