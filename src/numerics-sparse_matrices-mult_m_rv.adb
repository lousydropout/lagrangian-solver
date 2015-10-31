separate (Numerics.Sparse_Matrices)

function Mult_M_RV (Left  : in Sparse_Matrix;
		    Right : in Real_Vector) return Real_Vector is
   use Ada.Containers; 
   Vec   : Real_Vector := RV_Package.To_Vector (0.0, Count_Type (Left.N_Row));
   I     : Nat;
begin
   pragma Assert (Left.Format = CSC);

   for K in 1 .. Left.N_Col loop
      for J in Left.P (K) .. Left.P (K + 1) - 1 loop
	 I       := Left.I (J);
	 Vec (I) := Vec (I) + Left.X (J) * Right (K);
      end loop;
   end loop;
   return Vec;
end Mult_M_RV;
