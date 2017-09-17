with Ada.Text_IO;
use  Ada.Text_IO;
package Auto_Differentiation.Integrator is
   
   type Variable (N2 : Nat) is
      record
	 X : Real_Vector (1 .. N2);
	 T : Real;
      end record;
   
   type Control_Type is
      record
	 N   : Nat;
	 Dt  : Real := 1.0;
	 Eps : Real := 1.0e-10;
	 Err : Real := 1.0;
	 K   : Nat  := 9;
      end record;
   
   procedure FJ (Lagrangian : not null access 
		   function (X : Real_Vector; N : Nat) return AD_Type;
		 Var     : in     Variable;
		 Control : in Control_Type;
		 Q       : in     Real_Vector;
		 F       :    out Sparse_Vector;
		 J       :    out Sparse_Matrix);
   
   function Collocation (Lagrangian : not null access 
			   function (X : Real_Vector; N : Nat) return AD_Type;
			 Var        : in     Variable;
			 Control    : in out Control_Type) return Real_Vector
   with Pre => Is_Setup = True;
   
   function Bogack_Shampine (Hamiltonian : not null access 
			       function (X : Real_Vector; N : Nat) return AD_Type;
			     Var	 : in     Variable;
			     Control     : in out Control_Type)
			    return Real_Vector;
   
   procedure Update (Hamiltonian : not null access 
		    	       function (X : Real_Vector; N : Nat) return AD_Type;
		     Var         : in out Variable;
		     Control     : in out Control_Type);

   
   procedure Print_Data (Var : in Variable;
			 Hamiltonian : not null access 
			   function (X : Real_Vector; N : Nat) return AD_Type);
   
   procedure Print_XYZ (File : in File_Type;
			Var  : in Variable);
   
   
   procedure Print_Data_L (File	: in File_Type;
			   Var	: in Variable);
   
   
   procedure Setup (N : in Nat;
		    K : in Nat);
   
   function Is_Setup return Boolean;
   
private
   
   MatA, MatB, MatC, MatD : Sparse_Matrix;
   
end Auto_Differentiation.Integrator;
