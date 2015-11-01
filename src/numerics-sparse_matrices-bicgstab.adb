separate (Numerics.Sparse_Matrices)

function BiCGSTAB (A   : in     Sparse_Matrix;
		   B   : in     Real_Vector;
		   X0  : in     Real_Vector;
		   Err :    out Real;
		   Tol : in     Real	    := 1.0e-10) return Real_Vector is
   use RV_Package;
   R0, R, V, P, S, T, Res, X : Real_Vector;
   ρ, α, β, ω : Real := 1.0;
   Tmp : Real;
   I : Int := 1;
begin
   pragma Assert (A.N_Col = Nat (B.Length));
   X   := X0;
   R   := B - A * X;
   R0  := R;
   V   := To_Vector (0.0, B.Length);
   P   := V;
   Err := Norm2 (R);
   while Err > Tol loop
      Tmp := Dot_Product (R, R0);
      β   := (α / ω) * (Tmp / ρ);
      ρ   := Tmp;
      P   := R + β * (P - ω * V);
      V   := A * P;
      α   := ρ / Dot_Product (V, R0);
      S   := R - α * V;
      T   := A * S;
      ω   := Dot_Product (S, T) / Dot_Product (T, T);
      Res := α * P + ω * S;
      X   := X + Res;
      Err := Norm2 (Res);
      R   := S - ω * T;
   end loop;
   return X;
end BiCGSTAB;
