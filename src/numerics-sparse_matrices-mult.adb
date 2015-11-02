separate (Numerics.Sparse_Matrices)

function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix is
   use Ada.Containers;
   A  : Sparse_Matrix renames Left;
   B  : Sparse_Matrix renames Right;
   C  : Sparse_Matrix;
   X  : Real_Array (1 .. A.N_Row) := (others => 0.0);
   W  : Int_Array  (1 .. A.N_Row) := (others => 0);
   Nz : Pos := 1;
   N_Row : constant Count_Type := Count_Type (A.N_Row + 1);
   N_Res : constant Count_Type := Count_Type (A.N_Row * B.N_Col / 100);
   
begin
   
   C.Format := CSC; C.N_Row := A.N_Row; C.N_Col := B.N_Col;
   C.P.Reserve_Capacity (Count_Type (B.N_Col) + 1);
   C.I.Reserve_Capacity (N_Res);
   C.X.Reserve_Capacity (N_Res);
   
   for J in 1 .. B.N_Col loop
      if C.X.Capacity < C.X.Length + N_Row then
      	 C.I.Reserve_Capacity (C.X.Capacity + N_Row);
      	 C.X.Reserve_Capacity (C.X.Capacity + N_Row);
      end if;
      C.P.Append (Nz);
      
      for K in B.P (J) .. B.P (J + 1) - 1 loop
	 Scatter (A, B.I (K), B.X (K), W, X, J, C, Nz);
      end loop;
      
      for P in C.P (J) .. Nz - 1 loop
	 C.X.Append (X (C.I (P)));
      end loop;
   end loop;
   C.P.Append (Nz);
   
   -- Remove extra space hogged by sparse matrix C
   C.I.Reserve_Capacity (C.I.Length);
   C.X.Reserve_Capacity (C.X.Length);
   
   -- Need to sort entries
   C.Convert; C.Convert;
   
   return C;
end Mult;
