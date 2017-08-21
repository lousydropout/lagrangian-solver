with Numerics, Ada.Text_IO;
use  Numerics, Ada.Text_IO;

package Molecular_Dynamics is
   use Real_Functions, Real_IO, Int_IO;
   
   type BC_Vectype is array (Nat range <>) of Boolean;
   
   function Verlet (R	  : in Pos2D_Vector;
		    R_New : in Pos2D_Vector;
		    Is_BC : in BC_Vectype;
		    Dt	  : in Real) return Pos2D_Vector;

   function Calculate_Forces (R	    : in Pos2D_Vector;
			      Is_BC : in BC_Vectype) return Pos2D_Vector;
   
   function LJ_Force (From : in Pos2D;
		      To   : in Pos2D) return Pos2D;

   
   
   procedure Initialize_Lattice (N, M	 : in     Nat;
				 Vac	 : in     Nat;
				 Lattice :    out Pos2d_Vector;
				 Is_BC	 :    out BC_Vectype;
				 Sides	 : in     Boolean      := False);

   procedure Output (R	  : in Pos2D_Vector;
		     File : in File_Type);
   
   
   function Spring_Force (From : in Pos2D;
			  To   : in Pos2D) return Pos2D is (To - From);

      
end Molecular_Dynamics;
