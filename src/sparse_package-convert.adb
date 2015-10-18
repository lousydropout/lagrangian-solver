separate (Sparse_Package)

procedure Convert (Mat : in out Matrix) is
   X          : Real_Array (1 .. Nat (Mat.X.Length)) 
     := (others => 0.0);
   I          : Int_Array (1 .. Nat (Mat.I.Length))  
     := (others => 0);
   N          : constant Nat := Nat'Max (Mat.N_Col, Mat.N_Row);
   Row, Count : Int_Array (1 .. N + 1) := (others => 0);
   Index      : Int := 1;
   TRANSPOSE_EXCEPTION : exception;
begin
   case Mat.Format is
      when CSC => Mat.Format := CSR;
      when CSR => Mat.Format := CSC;
      when Triplet => raise TRANSPOSE_EXCEPTION;
   end case;
   
   for K of Mat.I loop
      Count (K) := Count (K) + 1;
   end loop;
   
   Row := Cumulative_Sum (Count); Count := Row;
   for K in 1 .. Nat (Mat.P.Length) - 1 loop
      for J in Mat.P (K) .. Mat.P (K + 1) - 1 loop
	 Index           := Row (Mat.I (J));
	 Row (Mat.I (J)) := Row (Mat.I (J)) + 1;
	 I (Index)       := K;
	 X (Index)       := Mat.X (J);
      end loop;
   end loop;
   
   case Mat.Format is
      when CSC => Mat.P := Vectorize (Count (1 .. Mat.N_Col + 1)); 
      when CSR => Mat.P := Vectorize (Count (1 .. Mat.N_Row + 1)); 
      when others => raise TRANSPOSE_EXCEPTION;
   end case;
   Mat.I := Vectorize (I); Mat.X := Vectorize (X);
end Convert;
