package body Numerics.Sparse_Matrices.CSparse is
   
   
   function LU_Decomposition (Mat : in Sparse_Matrix;
			      Tol : in Real   := 1.0e-20) return LU_Type is
      Sparse : Sparse_Ptr := To_Sparse (Mat);
      LU     : LU_Type;
   begin
      LU.Symbolic := CS_Sqr (Prob => Sparse);
      pragma Assert (LU.Symbolic /= null,
		     "ERROR: Problem with LU.Symbolic from CS_Sqr");
      LU.Numeric  := CS_LU (Sparse, LU.Symbolic, Creal (Tol));
      pragma Assert (LU.Numeric /= null,
		     "ERROR: Problem with LU.Numeric from CS_LU");
      LU.NCol     := Sparse.N;
      Sparse      := Free (Sparse);
      return LU;
   end LU_Decomposition;
   
   
   function Solve (LU  : in LU_Type;
		   B   : in Sparse_Vector;
		   Tol : in Real	  := 1.0e-20) return Sparse_Vector is
      Y : Real_Array := To_Array (B);
      X : Creal_Array (Cint (Y'First) .. Cint (Y'Last));
   begin
      for I in X'Range loop
	 X (I) := Creal (Y (Integer (I)));
      end loop;
      X := Solve (LU, X);
      for I in X'Range loop
	 Y (Integer (I)) := Real (X (I));
      end loop;
      pragma Assert (X'Length = B.NMax);
      return Sparse (Y, Tol => Tol);
   end Solve;

   
   function Solve (LU : in LU_Type;
		   B  : in Creal_Array) return Creal_Array is
      X : Creal_Ptrs.Pointer := Solve (LU, B);
      Y : Creal_Array (B'Range) with Convention => C, Address => X.all'Address;
   begin
      return Y;
   end Solve;
   
   function Solve (LU : in LU_Type;
		   B  : in Creal_Array) return Creal_Ptrs.Pointer is
      use C;
      Err : C.int;
      X : Creal_Ptrs.Pointer;
   begin
      X := Solve_CS (LU.NCol, LU.Symbolic, LU.Numeric, B, Err);
      pragma Assert (Err /= 1, "ERROR from Solve_CS");
      return X;
   end Solve;
      
   function Is_Valid (P	: in Creal_Ptrs.Pointer;
		      N	: in Cpos) return Boolean is
      X : Creal_Array (1 .. N) with Convention => C, Address => P.all'Address;
   begin
      return (for all Y of X => Y'Valid);
   end Is_Valid;
   
   
   function N_Col (LU : in LU_Type) return Pos is (Pos (LU.NCol));
   
   
   function To_Sparse (Mat : in Sparse_Matrix) return Sparse_Ptr is separate;
   
   
   
   procedure Free (LU : in out LU_Type) is
   begin
      LU.Symbolic := Free (LU.Symbolic);
      LU.Numeric  := Free (LU.Numeric);
      LU.NCol     := 0;
   end Free;
   
   function To_Array (X : in Real_Vector) return Creal_Array is
      Z : Real_Array := To_Array (X);
      Y : Creal_Array (Cint (Z'First) .. Cint (Z'Last));
   begin
      for I in Y'Range loop
	 Y (I) := Creal (Z (Integer (I)));
      end loop;
      return Y;
   end To_Array;
   
   function To_Array (X : in Int_Vector) return Cint_Array is
      Z : Int_Array := To_Array (X);
      Y : Cint_Array (Cint (Z'First) .. Cint (Z'Last));
   begin
      for I in Y'Range loop
	 Y (I) := Cint (Z (Integer (I)));
      end loop;
      return Y;
   end To_Array;
   
end Numerics.Sparse_Matrices.CSparse;

