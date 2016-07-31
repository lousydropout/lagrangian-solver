with Forward_AD, Forward_AD.Hamiltonian, Ada.Text_IO, Numerics;
use  Forward_AD, Forward_AD.Hamiltonian, Ada.Text_IO, Numerics;
--  with Forward_AD, Forward_AD.Hamiltonian, Ada.Text_IO, Numerics;
--  use  Forward_AD, Forward_AD.Hamiltonian, Ada.Text_IO, Numerics;

procedure AD_Test is
   
   A : AD_Type;
   
   Q : Real_Array := (1.0, 2.0, 3.0);
   V : Real_Array := Q;
   
   
   function Hamiltonian (Q : in out Real_Array) return Real is
   begin
      return 0.0;
   end Hamiltonian;
   
   
begin
   
   A := Var (30.0, 1, 2);
   --  B := A * A; -- Var (3.0, 2, 2);

   --  Put_Line ("------ A ----------------------");
   Print (A);
   --  Put_Line ("------ B ----------------------");
   --  Print (B);
   Put_Line ("------ A * B ----------------------");
   A := Func (Q, V, 0.0);
   Print (A);
   null;
end AD_Test;
