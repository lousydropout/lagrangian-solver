separate (Sparse_Package)

function BiCGSTAB (A	 : in Matrix;
		   B	 : in Real_Vector;
		   X0 : in Real_Vector) return Real_Vector is
   use RV_Package;
   R0, R, V, P, S, T, Res, X : Real_Vector;
   ρ, α, β, ω : Real := 1.0;
   Tmp, Err : Real;
begin
   pragma Assert (A.N_Col = Nat (B.Length));
   X   := X0;
   R   := B - A * X;
   R0  := R;
   V   := To_Vector (0.0, B.Length);
   P   := V;
   Err := Norm2 (R);
   while Err > 1.0e-20 loop
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
