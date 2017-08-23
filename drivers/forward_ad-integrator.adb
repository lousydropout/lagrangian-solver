with Numerics;
use  Numerics;

package body Forward_AD.Integrator is
   
   function Bogack_Shampine (Hamiltonian : not null access 
			       function (X : Real_Array; N : Nat) return AD_Type;
			     X		 : in     Real_Array;
			     N           : in     Nat;
			     T		 : in     Real;
			     Dt		 : in     Real;
			     Err	 :    out Real) 
			    return Real_Array is
      J : constant Sparse_Matrix := -Omega (N);
      K1, K2, K3, K4, Y, Z : Real_Array (X'Range);
   begin
      
      K1 := To_Array (J * Grad (Hamiltonian (X, N)));
      K2 := To_Array (J * Grad (Hamiltonian (X + (0.50 * Dt) * K1, N)));
      K3 := To_Array (J * Grad (Hamiltonian (X + (0.75 * Dt) * K2, N)));
      Y  := X + (Dt / 9.0) * (2.0 * K1 + 3.0 * K2 + 4.0 * K3);
      K4 := To_Array (J * Grad (Hamiltonian (Y, N)));
      Z  := X + (Dt / 24.0) * (7.0 * K1 + 6.0 * K2 + 8.0 * K3 + 3.0 * K4);

      Err := Norm (Z - Y);

      return (Z);
   end Bogack_Shampine;
   
   
   procedure Update (Hamiltonian : not null access 
		    	       function (X : Real_Array; N : Nat) return AD_Type;
		     X   : in out Real_Array;
		     N   : in     Nat;
		     T   : in out Real;
		     Dt  : in out Real;
		     Eps : in     Real  := 1.0e-10) is
      use Real_Functions;
      Y   : Real_Array (X'Range);
      Err : Real := 1.0;
   begin
      while Err > Eps loop
	 Y := Bogack_Shampine (Hamiltonian, X, N, T, Dt, Err);
	 if (Err <= Eps) then
	    X  := Y;
	    T  := T + Dt;
	 end if;
	 Dt := 0.8 * Dt * (Eps / Err) ** 0.3;
      end loop;
   end Update;
   
end Forward_AD.Integrator;
