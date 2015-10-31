separate (Numerics.Sparse_Matrices)

function Mult_M_RV (Left  : in Sparse_Matrix;
		    Right : in Real_Vector) return Real_Vector is
   use Ada.Containers, IV_Package; 
   V : Real_Array (1 .. Left.N_Row) := (others => 0.0);
   B : constant Real_Array  := To_Array (Right);
   I : Pos;
   L : Pos;
   R : Pos;
   C : Cursor;

begin
   pragma Assert (Left.Format = CSC);
   
   C := Left.P.First;
   R := Left.P (C);
   Next (C);
   for K in 1 .. Left.N_Col loop
      L := R;
      R := Left.P (C);
      Next (C);
      for J in L .. R - 1 loop
	 I     := Left.I (J);
	 V (I) := V (I) + Left.X (J) * B (K);
      end loop;
   end loop;
   return Vectorize (V);
end Mult_M_RV;
