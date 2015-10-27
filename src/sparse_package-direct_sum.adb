separate (Sparse_Package)

function Direct_Sum (Left, Right : in Sparse_Matrix) return Sparse_Matrix is
   Result : Sparse_Matrix;
   NRow : Pos := Pos (Left.N_Row);
   CounterI : Pos := Pos (Left.I.Length);
   CounterP : Pos := Pos (Left.P.Length);
   NP   : Pos := Left.P (CounterP);
   use Ada.Containers;
begin
   pragma Assert (Left.Format = CSC and Right.Format = CSC);
   Result.Format := CSC;
   Result.N_Col := Left.N_Col + Right.N_Col;
   Result.N_Row := Left.N_Row + Right.N_Row;
   Result.I.Set_Length (Left.I.Length + Right.I.Length);
   Result.P.Set_Length (Left.P.Length + Right.P.Length - 1);
   Result.X.Set_Length (Left.X.Length + Right.X.Length);

   -- left matrix
   for J in 1 .. Nat (Left.I.Length) loop
      Result.I (J) := Left.I (J);
      Result.X (J) := Left.X (J);
   end loop;

   for J in 1 .. Nat (Left.P.Length) loop
      Result.P (J) := Left.P (J);
   end loop;
   
   -- right matrix
   for J in 1 .. Pos (Right.I.Length) loop
      Result.I (CounterI + J) := NRow + Right.I (J);
      Result.X (CounterI + J) := Right.X (J);
   end loop;
   for J in 2 .. Pos (Right.P.Length) loop
      Result.P (CounterP + J - 1) := NP + Right.P (J) - 1;
   end loop;
   
   return Result;
end Direct_Sum;
