package body Numerics.Sparse_Matrices.CSparse is
   
   
   function LU_Decomposition (Mat : in Sparse_Matrix;
			      Tol : in Real   := 1.0e-12) return LU_Type is
      Sparse : Sparse_Ptr := To_Sparse (Mat);
      LU     : LU_Type;
   begin
      LU.Symbolic := CS_Sqr (Prob => Sparse);
      pragma Assert (LU.Symbolic /= null,
		     "ERROR: Problem with LU.Symbolic from CS_Sqr");
      LU.Numeric  := CS_LU (Sparse, LU.Symbolic, Tol);
      pragma Assert (LU.Numeric /= null,
		     "ERROR: Problem with LU.Numeric from CS_LU");
      LU.NCol     := Sparse.N;
      Sparse      := Free (Sparse);
      return LU;
   end LU_Decomposition;
   
   function Solve (LU  : in LU_Type;
		   B   : in Sparse_Vector;
		   Tol : in Real	  := 1.0e-20) return Sparse_Vector is
      X : Real_Array := Solve (LU, To_Array (B));
   begin
      pragma Assert (X'Length = B.NMax);
      return Sparse (X, Tol => Tol);
   end Solve;

   
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Array is
      X : Real_Ptrs.Pointer := Solve (LU, B);
      Y : Real_Array (B'Range) with Convention => C, Address => X.all'Address;
   begin
      return Y;
   end Solve;
   
   function Solve (LU : in LU_Type;
		   B  : in Real_Vector) return Real_Array is 
      (Solve (LU, To_Array (B)));
      
   function Solve (LU : in LU_Type;
		   B  : in Real_Vector) return Real_Vector is
      (Vectorize (Solve (LU, To_Array (B))));
      
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Ptrs.Pointer is
      use C;
      Err : C.int;
      X : Real_Ptrs.Pointer;
   begin
      X := Solve_CS (LU.NCol, LU.Symbolic, LU.Numeric, B, Err);
      pragma Assert (Err /= 1, "ERROR from Solve_CS");
      return X;
   end Solve;
      
   function Is_Valid (P	: in Real_Ptrs.Pointer;
		      N	: in Pos) return Boolean is
      X : Real_Array (1 .. N) with Convention => C, Address => P.all'Address;
   begin
      return (for all Y of X => Y'Valid);
   end Is_Valid;
   
   
   function N_Col (LU : in LU_Type) return Pos is (LU.NCol);
   
   
   function To_Sparse (Mat : in Sparse_Matrix) return Sparse_Ptr is separate;
   
   
   
   procedure Free (LU : in out LU_Type) is
   begin
      LU.Symbolic := Free (LU.Symbolic);
      LU.Numeric  := Free (LU.Numeric);
      LU.NCol     := 0;
   end Free;
end Numerics.Sparse_Matrices.CSparse;
