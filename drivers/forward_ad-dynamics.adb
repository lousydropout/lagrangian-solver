with Numerics, Ada.Text_IO, Forward_AD.Integrator;
use  Numerics, Ada.Text_IO, Forward_AD.Integrator;

procedure Forward_AD.Dynamics is
   use Real_Functions, Real_IO, Int_IO;
   function Hamiltonian (X : in Real_Array;
			 N : in Nat) return AD_Type;
   
   ---- Parameters -----
   N   : constant Nat  := 2;       -- # of particles
   Eps : constant Real := 1.0e-6;  -- error tolerance
   Dt  : Real := 1.0e2;            -- initial time-step
   -------------------------------				   
   
   --- Set up Hamiltonian -----
   function Hamiltonian (X : in Real_Array;
			 N : in Nat) return AD_Type is
      H : AD_Type   := Zero (X'Length);
      Q : AD_Vector := Var  (X (1     ..     N), 2 * N,     1);
      P : AD_Vector := Var  (X (N + 1 .. 2 * N), 2 * N, N + 1);
   begin
      for I in 1 .. N loop
	 H := H + 0.5 * P (I) ** 2;
      end loop;
      for I in 1 .. N - 1 loop
	 H := H + 0.5 * (Q (I + 1) - Q (I)) ** 2;
      end loop;
      return H;
   end Hamiltonian;
   -------------------------------
   
   
   -- Initial Conditions ----
   X   : Real_Array := 40.0 * Rand (2 * N);
   T   : Real       := 0.0;
   -------------------------------
begin
   
   Put (T); Put ("    "); Put (Dt); Put ("    ");
   Put (Val (Hamiltonian (X, N))); New_Line;
   
   Update (Hamiltonian'Access, X, N, T, Dt);
   
   Put (T); Put ("    "); Put (Dt); Put ("    ");
   Put (Val (Hamiltonian (X, N))); New_Line;


   
   null;
   
   Put_Line ("--------------------------------");
   for K in X'Range loop
      Put (X (K), Exp => 0, Aft => 3); New_Line;
   end loop;
   
end Forward_AD.Dynamics;
