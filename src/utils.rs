use ark_ff::PrimeField;
use ark_r1cs_std::{
    alloc::{AllocVar, AllocationMode},
    fields::nonnative::NonNativeFieldVar,
};
use ark_relations::r1cs::{Matrix as R1CSMatrix, Namespace, SynthesisError};
use core::{borrow::Borrow, marker::PhantomData};

pub struct SparseMatrix<F: PrimeField> {
    pub n_rows: usize,
    pub n_cols: usize,
    /// coeffs = R1CSMatrix = Vec<Vec<(F, usize)>>, which contains each row and the F is the value
    /// of the coefficient and the usize indicates the column position
    pub coeffs: R1CSMatrix<F>,
}

#[derive(Debug, Clone)]
pub struct SparseMatrixVar<F: PrimeField, CF: PrimeField> {
    _f: PhantomData<F>,
    _cf: PhantomData<CF>,
    pub n_rows: usize,
    pub n_cols: usize,
    // same format as the native SparseMatrix (which follows ark_relations::r1cs::Matrix format
    pub coeffs: Vec<Vec<(NonNativeFieldVar<F, CF>, usize)>>,
}

impl<F, CF> AllocVar<SparseMatrix<F>, CF> for SparseMatrixVar<F, CF>
where
    F: PrimeField,
    CF: PrimeField,
{
    fn new_variable<T: Borrow<SparseMatrix<F>>>(
        cs: impl Into<Namespace<CF>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        f().and_then(|val| {
            let cs = cs.into();

            let mut coeffs: Vec<Vec<(NonNativeFieldVar<F, CF>, usize)>> = Vec::new();
            for row in val.borrow().coeffs.iter() {
                let mut rowVar: Vec<(NonNativeFieldVar<F, CF>, usize)> = Vec::new();
                for &(value, col_i) in row.iter() {
                    let coeffVar =
                        NonNativeFieldVar::<F, CF>::new_variable(cs.clone(), || Ok(value), mode)?;
                    rowVar.push((coeffVar, col_i));
                }
                coeffs.push(rowVar);
            }

            Ok(Self {
                _f: PhantomData,
                _cf: PhantomData,
                n_rows: val.borrow().n_rows,
                n_cols: val.borrow().n_cols,
                coeffs,
            })
        })
    }
}

pub fn mat_vec_mul_sparse<F: PrimeField>(M: &SparseMatrix<F>, z: &[F]) -> Vec<F> {
    assert_eq!(M.n_cols, z.len());
    let mut res = vec![F::zero(); M.n_rows];
    for (row_i, row) in M.coeffs.iter().enumerate() {
        for &(value, col_i) in row.iter() {
            res[row_i] += value * z[col_i];
        }
    }
    res
}
pub fn dense_matrix_to_sparse<F: PrimeField>(m: Vec<Vec<F>>) -> SparseMatrix<F> {
    let mut r = SparseMatrix::<F> {
        n_rows: m.len(),
        n_cols: m[0].len(),
        coeffs: Vec::new(),
    };
    for m_row in m.iter() {
        let mut row: Vec<(F, usize)> = Vec::new();
        for (col_i, value) in m_row.iter().enumerate() {
            if !value.is_zero() {
                row.push((*value, col_i));
            }
        }
        r.coeffs.push(row);
    }
    r
}

// just some helpers to define matrices and vectors by hand
pub fn to_F_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> SparseMatrix<F> {
    dense_matrix_to_sparse(to_F_dense_matrix(M))
}
pub fn to_F_dense_matrix<F: PrimeField>(M: Vec<Vec<usize>>) -> Vec<Vec<F>> {
    M.iter()
        .map(|m| m.iter().map(|r| F::from(*r as u64)).collect())
        .collect()
}
pub fn to_F_vec<F: PrimeField>(z: Vec<usize>) -> Vec<F> {
    z.iter().map(|c| F::from(*c as u64)).collect()
}
