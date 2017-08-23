package Forward_AD.Integrator is
   
   function Bogack_Shampine (Hamiltonian : not null access 
			       function (X : Real_Array; N : Nat) return AD_Type;
			     X		 : in     Real_Array;
			     N           : in     Nat;
			     T		 : in     Real;
			     Dt		 : in     Real;
			     Err	 :    out Real)
			    return Real_Array;
   
   procedure Update (Hamiltonian : not null access 
		    	       function (X : Real_Array; N : Nat) return AD_Type;
		     X   : in out Real_Array;
		     N   : in     Nat;
		     T   : in out Real;
		     Dt  : in out Real;
		     Eps : in     Real  := 1.0e-10);


   
end Forward_AD.Integrator;
