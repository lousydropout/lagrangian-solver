separate (Numerics.Sparse_Matrices)

function Direct_Sum (A, B : in Sparse_Matrix) return Sparse_Matrix is
   C    : Sparse_Matrix;
   NRow : constant Pos := A.N_Row;
   NP   : constant Pos := A.P (A.N_Col + 1);
   use Ada.Containers;
begin
   pragma Assert (A.Format = CSC and B.Format = CSC);
   C.Format := CSC;
   C.N_Col  := A.N_Col + B.N_Col;
   C.N_Row  := A.N_Row + B.N_Row;
   C.X.Reserve_Capacity (A.X.Length + B.X.Length);
   C.I.Reserve_Capacity (A.I.Length + B.I.Length);
   C.P.Reserve_Capacity (A.P.Length + B.P.Length - 1);

   -- left matrix
   for J in 1 .. Pos (A.I.Length) loop
      C.X.Append (A.X (J));
   end loop;
   for J in 1 .. Pos (A.I.Length) loop
      C.I.Append (A.I (J));
   end loop;
   for J in 1 .. Pos (A.P.Length) loop
      C.P.Append (A.P (J));
   end loop;
   
   -- right matrix
   for J in 1 .. Pos (B.I.Length) loop
      C.X.Append (B.X (J));
   end loop;
   for J in 1 .. Pos (B.I.Length) loop
      C.I.Append (NRow + B.I (J));
   end loop;
   for J in 2 .. Pos (B.P.Length) loop
      C.P.Append (NP + B.P (J) - 1);
   end loop;
   
   return C;
end Direct_Sum;

