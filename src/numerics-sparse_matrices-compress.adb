separate (Numerics.Sparse_Matrices)

procedure Compress (Mat : in out Sparse_Matrix) is
   X          : Real_Array (1 .. Int (Mat.X.Length)) := (others => 0.0);
   I          : Int_Array  (1 .. Int (Mat.X.Length)) := (others => 0);
   Col, Count : Int_Array  (1 .. Mat.N_Col + 1)      := (others => 0);
   Index      : Nat                                  := 1;
   P          : Int_Vector renames Mat.P;
begin
   Mat.Format := CSC;
   
   for K of P loop
      Count (K) := Count (K) + 1;
   end loop;
   
   Col := Cumulative_Sum (Count); Count := Col;

   for K in 1 .. Nat (Mat.X.Length) loop
      Index       := Col (P (K));
      Col (P (K)) := Col (P (K)) + 1;
      I (Index)   := Mat.I (K); 
      X (Index)   := Mat.X (K);
   end loop;
   
   Set (Mat.I, I);
   Set (Mat.X, X);
   Set (Mat.P, Count);
   
   Mat.Convert; Mat.Convert;
   Mat.Remove_Duplicates;
end Compress;
