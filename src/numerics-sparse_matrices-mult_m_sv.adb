separate (Numerics.Sparse_Matrices)

function Mult_M_SV (A : in Sparse_Matrix;
		    X : in Sparse_Vector) return Sparse_Vector is 
   I : constant Int_Array  := To_Array (X.I);
   B : constant Real_Vector := To_Array (X.X);
   V : Real_Vector (1 .. A.N_Row) := (others => 0.0);
   N : Pos;
   
begin
   pragma Assert (A.Format = CSC);
   pragma Assert (I'Length = Pos (X.I.Length));
   
   for K in I'Range loop
      for J in A.P (I (K)) .. A.P (I (K) + 1) - 1 loop
   	 N     := A.I (J);
   	 V (N) := V (N) + A.X (J) * B (K);
      end loop;
   end loop;
   
   return Sparse (V, Tol => 1.0e-15);
end Mult_M_SV;
