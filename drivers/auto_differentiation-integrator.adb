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

      Level := Hessian; -- Turn back on calculation of Hessians
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
	 Dt := 0.8 * Dt * (Eps / Err) ** 0.3;
      end loop;
   end Update;
   
end Auto_Differentiation.Integrator;
