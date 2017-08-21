with Numerics;
use  Numerics;

package body Forward_AD.Integrator is
   
   function Bogack_Shampine (X		 : in     Real_Array;
			     T		 : in     Real;
			     Dt		 : in     Real;
			     Err	 :    out Real) return Real_Array is
      J : constant Sparse_Matrix := -Omega (N);
      K1, K2, K3, K4, Y, Z : Real_Array (X'Range);
   begin
      
      K1 := To_Array (J * Grad (Hamiltonian (X)));
      K2 := To_Array (J * Grad (Hamiltonian (X + (0.50 * Dt) * K1)));
      K3 := To_Array (J * Grad (Hamiltonian (X + (0.75 * Dt) * K2)));
      Y  := X + (Dt / 9.0) * (2.0 * K1 + 3.0 * K2 + 4.0 * K3);
      K4 := To_Array (J * Grad (Hamiltonian (Y)));
      Z  := X + (Dt / 24.0) * (7.0 * K1 + 6.0 * K2 + 8.0 * K3 + 3.0 * K4);

      Err := Norm (Z - Y);

      return (Z);
   end Bogack_Shampine;
   
   
   procedure Update (X		 : in out Real_Array;
		     T		 : in out Real;
		     Dt		 : in out Real) is
      use Real_Functions;
      Y   : Real_Array (X'Range);
      Err : Real := 2.0 * Eps;
   begin
      while Err > Eps loop
	 Y := Bogack_Shampine (X, T, Dt, Err);
	 if (Err <= Eps) then
	    X  := Y;
	    T  := T + Dt;
	 end if;
	 Dt := 0.8 * Dt * (Eps / Err) ** 0.3;
      end loop;
   end Update;
   
end Forward_AD.Integrator;
