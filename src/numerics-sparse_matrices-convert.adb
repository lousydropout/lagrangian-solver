separate (Numerics.Sparse_Matrices)

procedure Convert (Mat : in out Sparse_Matrix) is
   N     : constant Pos := Pos'Max (Mat.N_Col, Mat.N_Row);
   Row   : Int_Array  (1 .. N + 1)   := (others => 0);
   Count : Int_Array  (1 .. N + 1)   := (others => 0);
   Nmax  : constant Pos := Pos (Mat.X.Length);
   X     : Real_Vector (1 .. Nmax)    := (others => 0.0);
   I     : Int_Array  (1 .. Nmax)    := (others => 0);
   Index : Nat                       := 1;
   Tmp   : Nat                       := 1;
   
   Transpose_Exception : exception;
begin
   case Mat.Format is
      when CSC => Mat.Format := CSR;
      when CSR => Mat.Format := CSC;
      when Triplet => raise Transpose_Exception;
   end case;
   
   for K of Mat.I loop Count (K) := Count (K) + 1; end loop;
   
   Cumulative_Sum (Count); Row := Count;

   for K in 1 .. Nat (Mat.P.Length) - 1 loop
      for J in Mat.P (K) .. Mat.P (K + 1) - 1 loop
	 Tmp       := Mat.I (J);
	 Index     := Row (Tmp);
	 Row (Tmp) := Row (Tmp) + 1;
	 I (Index) := K;
	 X (Index) := Mat.X (J);
      end loop;
   end loop;
   
   case Mat.Format is
      when CSC => Set (Mat.P, Count (1 .. Mat.N_Col + 1)); 
      when CSR => Set (Mat.P, Count (1 .. Mat.N_Row + 1)); 
      when others => raise Transpose_Exception;
   end case;
   Set (Mat.I, I); 
   Set (Mat.X, X);
end Convert;
