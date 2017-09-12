with Numerics;
use  Numerics;

package body Auto_Differentiation.Integrator is
   
   function Bogack_Shampine (Hamiltonian : not null access 
			       function (X : Real_Array; N : Nat) return AD_Type;
			     Var	 : in     Variable;
			     Control     : in out Control_Type)
			    return Real_Array is
      X   : Real_Array renames Var.X;
      N   : Nat  renames Control.N;
      Dt  : Real renames Control.Dt;
      Err : Real renames Control.Err;
      J   : constant Sparse_Matrix := -Omega (N);
      K1, K2, K3, K4, Y, Z : Real_Array (X'Range);
      Old : constant Evaluation_Level := Level;
   begin
      pragma Assert (2 * N = Var.N2);
      -- Turn off the calculation of Hessians (not used for explicit schemes):
      Level := Gradient; 
      
      K1 := To_Array (J * Grad (Hamiltonian (X, N)));
      K2 := To_Array (J * Grad (Hamiltonian (X + (0.50 * Dt) * K1, N)));
      K3 := To_Array (J * Grad (Hamiltonian (X + (0.75 * Dt) * K2, N)));
      Y  := X + (Dt / 9.0) * (2.0 * K1 + 3.0 * K2 + 4.0 * K3);
      K4 := To_Array (J * Grad (Hamiltonian (Y, N)));
      Z  := X + (Dt / 24.0) * (7.0 * K1 + 6.0 * K2 + 8.0 * K3 + 3.0 * K4);

      Err := Norm (Z - Y);

      Level := Old; -- return to previous evaluation level
      return (Z);
   end Bogack_Shampine;
   
   
   procedure Update (Hamiltonian : not null access 
		    	       function (X : Real_Array; N : Nat) return AD_Type;
		     Var         : in out Variable;
		     Control     : in out Control_Type) is
      use Real_Functions;
      X   : Real_Array renames Var.X;
      T   : Real renames Var.T;
      N   : Nat  renames Control.N;
      Dt  : Real renames Control.Dt;
      Err : Real renames Control.Err;
      Eps : Real renames Control.Eps;
      Y   : Real_Array (X'Range);
   begin
      pragma Assert (2 * N = Var.N2);
      
      Err := 1.0;
      while Err > Eps loop
	 Y := Bogack_Shampine (Hamiltonian, Var, Control);
	 if (Err <= Eps) then
	    X  := Y;
	    T  := T + Dt;
	 end if;
	 Dt := 0.8 * Dt * (Eps / (Err + 1.0e-20)) ** 0.3;
      end loop;
   end Update;
   
   
   procedure Print_XYZ (File : in File_Type;
			Var  : in Variable) is
      use Real_Functions, Real_IO;
      X :  Real_Array renames Var.X;
      X1, Y1, X2, Y2 : Real;
      R : constant Real := 10.0;
   begin
      X1 := -R * Sin (X (1));
      Y1 :=  R * Cos (X (1));
      X2 := X1 - R * Sin (2.0 * X (1) + X (2));
      Y2 := Y1 + R * Cos (2.0 * X (1) + X (2));
      -- print header
      Put_Line (File, "3");
      Put (File, "Properties=pos:R:2   Time=");
      Put (File, Var.T, Aft => 5, Exp => 0);
      New_Line(File);
      -- position of ball 1
      Put_Line (File, "0.0     0.0     5.0");
      -- position of ball 2
      Put (File => File, Item => X1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
      -- position of ball 3
      Put (File => File, Item => X2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
   end Print_XYZ;
   
   --- print data ------
   procedure Print_Data (Var : in Variable;
			 Hamiltonian : not null access 
			   function (X : Real_Array; N : Nat) return AD_Type) is
      use Real_Functions, Real_IO;
      T : Real renames Var.X (1);
      S : Real renames Var.X (2);
      X, Y : Real_Array (1 .. 2);
   begin
      X (1) := -Sin (T);
      Y (1) :=  Cos (T);
      X (2) := X (1) - Sin (2.0 * T + S);
      Y (2) := Y (1) + Cos (2.0 * T + S);
      
      ---------------------------------------
      Put (Var.T, Aft => 6, Exp => 0); -- print time
      for I in 1 .. 2 loop
	 Put (",  "); Put (X (I), Aft => 4, Exp => 0);
	 Put (",  "); Put (Y (I), Aft => 4, Exp => 0); -- print positions
      end loop;
      Put (",  "); 
      -- print total energy
      Put (Val (Hamiltonian (Var.X, 2)), Aft => 10, Exp => 0); New_Line;
   end Print_Data;
   
   
end Auto_Differentiation.Integrator;
