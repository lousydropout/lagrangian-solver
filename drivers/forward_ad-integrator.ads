generic
   
   N   : Nat;
   Eps : Real;
   with function Hamiltonian (X : in Real_Array) return AD_Type;
   
package Forward_AD.Integrator is

   function Bogack_Shampine (X		 : in     Real_Array;
			     T		 : in     Real;
			     Dt		 : in     Real;
			     Err	 :    out Real) return Real_Array;
   
   procedure Update (X		 : in out Real_Array;
		     T		 : in out Real;
		     Dt		 : in out Real);


   
end Forward_AD.Integrator;
