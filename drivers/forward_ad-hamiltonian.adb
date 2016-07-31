package body Forward_AD.Hamiltonian is
   
   function Func (Q, V : in Real_Array; T : in Real) return AD_Type is
      N : constant Nat := Q'Length;
   begin
      pragma Assert (Q'Length = V'Length);
      return Var (1.0, 1, 1);
   end Func;

   
end Forward_AD.Hamiltonian;
