with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

package Molecular_Dynamics is
   use Real_Functions, Real_IO, Int_IO;
   
   
   function Verlet (R	  : in Pos2D_Vector;
		    R_New : in Pos2D_Vector;
		    Dt	  : in Real) return Pos2D_Vector;

   
   function Correct_Forces (F : in Pos2D_Vector;
			    M : in Nat) return Pos2D_Vector;
   
   function Calculate_Forces (R : in Pos2D_Vector) return Pos2D_Vector;
   
   function LJ_Force (From : in Pos2D;
		      To   : in Pos2D) return Pos2D;

   
   
   function Initialize_Lattice (N, M : in Nat;
				Vac  : in Nat) return Pos2D_Vector;

   procedure Output (R	  : in Pos2D_Vector;
		     File : in File_Type);
   
   
   function Spring_Force (From : in Pos2D;
			  To   : in Pos2D) return Pos2D is (To - From);

      
end Molecular_Dynamics;
