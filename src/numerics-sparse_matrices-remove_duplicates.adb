separate (Numerics.Sparse_Matrices)

procedure Remove_Duplicates (Mat : in out Sparse_Matrix) is
   N, Iter : Pos := 0;
   J : Int_Array  (1 .. Nat (Mat.P.Length)) := (others => 0);
   I : Int_Array  (1 .. Nat (Mat.I.Length));
   X : Real_Vector (1 .. Nat (Mat.X.Length));
begin
   for K in 1 .. Nat (Mat.P.Length) - 1 loop
      Iter := 0;
      for L in Mat.P (K) .. Mat.P (K + 1) - 1 loop
	 if Iter /= Mat.I (L) then
	    N     := N + 1;
	    Iter  := Mat.I (L);
	    I (N) := Iter; 
	    J (K) := J (K) + 1;
	    X (N) := Mat.X (L);
	 else
	    X (N) := X (N) + Mat.X (L);
	 end if;
      end loop;
   end loop;
   Cumulative_Sum (J);
   Set (Mat.I, I (1 .. N));
   Set (Mat.X, X (1 .. N));
   Set (Mat.P, J);
end Remove_Duplicates;
