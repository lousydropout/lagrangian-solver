with Numerics, Ada.Text_IO, Auto_Differentiation.Integrator;
use  Numerics, Ada.Text_IO, Auto_Differentiation.Integrator;

procedure Auto_Differentiation.Dynamics is
   use Real_IO, Int_IO;
   --  Set Up Parameters -----------------
   N       : constant Nat := 2;
   Control : Control_Type (N => N);
   -------------------------------
   
   --- Set up Hamiltonian -----
   function Hamiltonian (X : in Real_Array; N : in Nat) return AD_Type is
      H : AD_Type   := Zero (X'Length);
      Q : AD_Vector := Var  (X (1     ..     N), 2 * N,     1);
      P : AD_Vector := Var  (X (N + 1 .. 2 * N), 2 * N, N + 1);
      Z : AD_Type;
   begin
      for I in 1 .. N loop
	 H := H + 0.5 * P (I) ** 2;
      end loop;
      for I in 1 .. N - 1 loop
	 Z := Q (I + 1) - Q (I);
	 H := H + 0.5 * (Z ** 2);
      end loop;
      return H;
   end Hamiltonian;
   -------------------------------
   
   -- Initial Conditions ----
   Var : Variable := (N2 => 2 * N,
		      X => 40.0 * Rand (2 * N), 
		      T => 0.0);
   X : Real_Array renames Var.X;
   T : Real       renames Var.T;
   -------------------------------
begin
   Level := Value;
   
   --  Put (T, Fore => 3, Exp => 0, Aft => 3); Put ("    "); 
   --  Put (Control.Dt); Put ("    ");
   --  Put (Val (Hamiltonian (X, N))); New_Line;
   
   for Iter in 1 .. 4 loop
      for Iter2 in 1 .. 1000 loop
	 Update (Hamiltonian'Access, Var, Control);
      end loop;

      --  Put (T, Fore => 3, Exp => 0, Aft => 3); Put ("    "); 
      --  Put (Control.Dt); Put ("    ");
      --  Put (Val (Hamiltonian (X, N))); New_Line;
   end loop;
   
   null;
   
   --  Put_Line ("--------------------------------");
   --  for K in X'Range loop
   --     Put (X (K), Exp => 0, Aft => 3); New_Line;
   --  end loop;
   
end Auto_Differentiation.Dynamics;
