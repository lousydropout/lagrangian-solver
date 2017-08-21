with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

procedure Forward_AD.Verlet is
   use Real_Functions, Real_IO, Int_IO;
   
   --  function Potential (X : in Real_Array) return AD_Type;
   --  function Potential (X : in Real_Array) return AD_Type is
   --     Q : AD_Vector := Var (X);
   --     U : AD_Type   := Zero (X'Length);
   --  begin
   --     for I in Q'Range loop
   --  	 U := U + Q (I) ** 2;
   --     end loop;
   --     return U;
   --  end Potential;
   
   
   
   function Hamiltonian (X : in Real_Array;
			 N : in Nat) return AD_Type
     with Pre => X'Length = 2 * N;
   function Hamiltonian (X : in Real_Array;
			 N : in Nat) return AD_Type is
      Y : AD_Vector := Var (X);
      H : AD_Type   := Zero (X'Length);
      Q : AD_Vector := Y (1 .. N);
      P : AD_Vector (1 .. N);
      K : constant Real := 1.0;
   begin
      for I in 1 .. N loop
	 P (I) := Y (N + I);
      end loop;
      
      for I in 1 .. N loop
	 H := H + 0.5 * P (I) ** 2;
      end loop;
      for I in 1 .. N - 1 loop
	 H := H + 0.5 * K * (Q (I + 1) - Q (I)) ** 2;
      end loop;
      
      return H;
   end Hamiltonian;
   
   X : Real_Array := 4.0 * Rand (4);
   H : AD_Type := Hamiltonian (X, 2);
   --  F : Sparse_Vector := -H.Grad;
   F : Real_Array := -To_Array (Omega (2) * Grad (H));
   G, D : Real_Array (1 .. 4);
   --  F : Real_Array := -Grad (Potential (X));

begin

   G (1) := X (3);
   G (2) := X (4);
   G (3) := X (2) - X (1);
   G (4) := -G (3);
   
   D := F - G;
   null;
   Print (H);
   Put_Line ("--------------------------------");

   Put_Line ("--------------------------------");
   for K in F'Range loop
      Put (F (K), Exp => 0, Aft => 3); 
      Put ("    ");
      Put (G (K), Exp => 0, Aft => 3); 
      Put ("    ");
      Put (D (K), Exp => 0, Aft => 3); 
      New_Line;
   end loop;
   
end Forward_AD.Verlet;
