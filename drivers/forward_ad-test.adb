with Ada.Text_IO, Numerics;
use  Ada.Text_IO, Numerics;

procedure Forward_AD.Test is
   
   function Hamiltonian (Q, V : in Real_Array; T : in Real) return AD_Type is
      N : constant Nat := Q'Length;
   begin
      pragma Assert (Q'Length = V'Length);
      return Var (1.0, 1, 1);
   end Hamiltonian;
   
   
   A : AD_Type;
   
   Q : Real_Array := (1.0, 2.0, 3.0);
   V : Real_Array := Q;
   
begin
   A := Hamiltonian (Q, V, 0.0);
   Print (A);
   null;
end Forward_AD.Test;
