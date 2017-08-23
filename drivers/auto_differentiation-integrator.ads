package Auto_Differentiation.Integrator is
   
   type Variable (N2 : Nat) is
      record
	 X : Real_Array (1 .. N2);
	 T : Real;
      end record;
   
   type Control_Type (N : Nat) is
      record
	 Dt  : Real := 1.0;
	 Eps : Real := 1.0e-10;
	 Err : Real := 1.0;
      end record;
   
   function Bogack_Shampine (Hamiltonian : not null access 
			       function (X : Real_Array; N : Nat) return AD_Type;
			     Var	 : in     Variable;
			     Control     : in out Control_Type)
			    return Real_Array;
   
   procedure Update (Hamiltonian : not null access 
		    	       function (X : Real_Array; N : Nat) return AD_Type;
		     Var         : in out Variable;
		     Control     : in out Control_Type);

   
end Auto_Differentiation.Integrator;
