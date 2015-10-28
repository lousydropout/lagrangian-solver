separate (Numerics.Sparse_Matrices)

function Mult_M_RV (Left  : in Sparse_Matrix;
		    Right : in Real_Vector) return Real_Vector is
   Vec : Real_Vector;
   Mat : Sparse_Matrix;
   Tmp : Real;
   Index : Nat;
begin
   pragma Assert (Left.Format = CSC or Left.Format = CSR);

   Vec.Set_Length (Right.Length);

   if Left.Format = CSC then
      Mat := Convert (Left);
   else
      Mat := Left;
   end if;
   
   for K in 1 .. Left.N_Row loop
      Tmp := 0.0;
      for J in Mat.P (K) .. Mat.P (K + 1) - 1 loop
	 Index := Mat.I (J);
	 Tmp   := Tmp + Mat.X (J) * Right (Index);
      end loop;
      Vec (K) := Tmp;
   end loop;
   return Vec;
end Mult_M_RV;
