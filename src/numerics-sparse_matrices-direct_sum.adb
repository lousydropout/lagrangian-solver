separate (Numerics.Sparse_Matrices)

function Direct_Sum (A, B : in Sparse_Matrix) return Sparse_Matrix is
   C : Sparse_Matrix;
   Tmp  : Int_Vector;
   NRow : constant Pos := A.N_Row;
   NP   : constant Pos := A.P.Last_Element - 1;
   use Ada.Containers;
begin
   pragma Assert (A.Format = CSC and B.Format = CSC);
   C.Format := CSC;
   C.N_Col := A.N_Col + B.N_Col;
   C.N_Row := A.N_Row + B.N_Row;
   C.X.Reserve_Capacity (A.X.Length + B.X.Length);
   C.I.Reserve_Capacity (A.I.Length + B.I.Length);
   C.P.Reserve_Capacity (A.P.Length + B.P.Length - 1);
   Tmp.Reserve_Capacity (B.P.Length);

   for X of A.X loop C.X.Append (X); end loop;
   for X of B.X loop C.X.Append (X); end loop;
   
   for I of A.I loop C.I.Append (I); end loop;
   for I of B.I loop C.I.Append (NRow + I); end loop;
   
   for P of A.P loop C.P.Append (P); end loop;
   for P of reverse B.P loop Tmp.Append (P); end loop;
   Tmp.Delete_Last;
   for P of Tmp loop C.P.Append (NP + P); end loop;
   
   return C;
end Direct_Sum;


--  function Direct_Sum (Left, Right : in Sparse_Matrix) return Sparse_Matrix is
--     Result : Sparse_Matrix;
--     NRow : Pos := Pos (Left.N_Row);
--     CounterI : Pos := Pos (Left.I.Length);
--     CounterP : Pos := Pos (Left.P.Length);
--     NP   : Pos := Left.P (CounterP);
--     use Ada.Containers;
--  begin
--     pragma Assert (Left.Format = CSC and Right.Format = CSC);
--     Result.Format := CSC;
--     Result.N_Col := Left.N_Col + Right.N_Col;
--     Result.N_Row := Left.N_Row + Right.N_Row;
--     Result.I.Set_Length (Left.I.Length + Right.I.Length);
--     Result.P.Set_Length (Left.P.Length + Right.P.Length - 1);
--     Result.X.Set_Length (Left.X.Length + Right.X.Length);

--     -- left matrix
--     for J in 1 .. Nat (Left.I.Length) loop
--        Result.I (J) := Left.I (J);
--        Result.X (J) := Left.X (J);
--     end loop;

--     for J in 1 .. Nat (Left.P.Length) loop
--        Result.P (J) := Left.P (J);
--     end loop;
   
--     -- right matrix
--     for J in 1 .. Pos (Right.I.Length) loop
--        Result.I (CounterI + J) := NRow + Right.I (J);
--        Result.X (CounterI + J) := Right.X (J);
--     end loop;
--     for J in 2 .. Pos (Right.P.Length) loop
--        Result.P (CounterP + J - 1) := NP + Right.P (J) - 1;
--     end loop;
   
--     return Result;
--  end Direct_Sum;

