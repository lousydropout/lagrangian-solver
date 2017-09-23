with Numerics, Numerics.Sparse_Matrices, Chebyshev, Ada.Text_IO;
use  Numerics, Numerics.Sparse_Matrices, Chebyshev, Ada.Text_IO;
generic
   K : in Nat;
   N_Constraints : in Pos := 0;
package Dense_AD.Integrator is
   
   N : constant Nat := (Num - N_Constraints) / 2;

   type Dense_Or_Sparse is (Dense, Sparse);
   type Create_Or_Append is (Create, Append);
   
   type Variable is record
      X : Vector;
      T : Real;
   end record;
   
   type Control_Type is record
      Dto     : Real := 1.0;
      Dt      : Real := 1.0;
      Dtn     : Real := 1.0;
      Tol     : Real := 1.0e-10;
      Err     : Real := 1.0;
      Started : Boolean := False;
   end record;
   
   type Array_Of_Vectors is array (1 .. Num) of Real_Vector (1 .. K);
   
   function New_Control_Type (Dt      : in Real	   := 1.0;
			      Tol     : in Real	   := 1.0e-10;
			      Started : in Boolean := False) return Control_Type;

   
   function Chebyshev_Transform (Y : in Real_Vector) return Array_Of_Vectors
     with Pre => Y'First = 1 and Y'Length = Num * K;
   
   
   function Interpolate (A : in Array_Of_Vectors;
			 T : in Real;
			 L : in Real;
			 R : in Real) return Vector;
   
   procedure Update (Var : in out Variable;
		     Y	 : in     Real_Vector;
		     Dt	 : in     Real)
     with Pre => Y'First = 1 and Y'Length = Num * K;
   
   function Update (Lagrangian : not null access 
		      function (T : real; X : Vector) return AD_Type;
		    Var        : in     Variable;
		    Control    : in out Control_Type;
		    Density    : in Dense_Or_Sparse) return Real_Vector;
      
   function Update (Lagrangian : not null access 
		      function (T : real; X : Vector) return AD_Type;
		    Var        : in     Variable;
		    Control    : in out Control_Type) return Real_Vector;
      
   procedure Print_Lagrangian (File	  : in     File_Type;
			       Var	  : in     Variable;
			       Lagrangian : not null access
				 function (T : real; X : Vector) return AD_Type;
			       Fore : in Field := 3;
			       Aft  : in Field := 5;
			       Exp  : in Field := 3);
   
   procedure Print_Lagrangian (Var	  : in     Variable;
			       Lagrangian : not null access
				 function (T : real; X : Vector) return AD_Type;
			       Fore : in Field := 3;
			       Aft  : in Field := 5;
			       Exp  : in Field := 3);
   
   procedure Print_XYZ (File : in out File_Type;
			Var  : in     Variable);
   
   procedure Print_XYZ (File : in out File_Type;
			Var  : in     Variable;
			Name : in     String;
			Mode : in     Create_Or_Append := Append);
   
private

   Collocation_Density : Dense_Or_Sparse := Sparse;
   
   function Split (Y : in Real_Vector) return Array_Of_Vectors
     with Pre => Y'First = 1 and Y'Length = Num * K;
   
   procedure Iterate (Lagrangian : not null access 
			function (T : real; X : Vector) return AD_Type;
		      Y          : in out Real_Vector;
		      Var        : in     Variable;
		      Control    : in out Control_Type);
   
   procedure Colloc (Lagrangian : not null access 
			  function (T : real; X : Vector) return AD_Type;
			Q          : in out Real_Vector;
			Var        : in     Variable;
			Control    : in out Control_Type);
   
   procedure Collocation (Lagrangian : not null access 
			    function (T : real; X : Vector) return AD_Type;
			  Q          : in out Real_Vector;
			  Var        : in     Variable;
			  Control    : in out Control_Type);
   procedure FJ (Lagrangian : not null access 
		    function (T : real; X : Vector) return AD_Type;
		  Var     : in     Variable;
		  Control : in     Control_Type;
		  Q       : in     Real_Vector;
		  F       :    out Real_Vector;
		  J       :    out Real_Matrix);
   
   procedure Sp_Collocation (Lagrangian : not null access 
			       function (T : real; X : Vector) return AD_Type;
			     Q          : in out Real_Vector;
			     Var        : in     Variable;
			     Control    : in out Control_Type);
   procedure Sp_FJ (Lagrangian : not null access 
		      function (T : real; X : Vector) return AD_Type;
		    Var     : in     Variable;
		    Control : in     Control_Type;
		    Q       : in     Real_Vector;
		    F       :    out Sparse_Vector;
		    J       :    out Sparse_Matrix);
   
   type Array_Of_Sparse_Matrix is array (1 .. 4) of Sparse_Matrix;
   
   function Setup return Array_Of_Sparse_Matrix;
   
   Grid   : constant Real_Vector := Chebyshev_Gauss_Lobatto (K, 0.0, 1.0);
   Der    : constant Real_Matrix := Derivative_Matrix (K, 0.0, 1.0);
   NK     : constant Nat := Num * K;
   
   EyeN         : constant Sparse_Matrix := Eye (N);
   EyeK         : constant Sparse_Matrix := Eye (K);
   Eye2N        : constant Sparse_Matrix := Eye (2 * N);
   EyeNc        : constant Sparse_Matrix := Eye (N_Constraints);
   D            : constant Sparse_Matrix := Sparse (Der);
   ZeroN        : constant Sparse_Matrix := Zero (N);
   ZeroNc       : constant Sparse_Matrix := Zero (N_Constraints);
   Zero2N       : constant Sparse_Matrix := Zero (2 * N);
   
   Sp    : constant Array_Of_Sparse_Matrix := Setup;
   Mat_A : constant Real_Matrix := Dense (Sp (1));
   Mat_B : constant Real_Matrix := Dense (Sp (2));
   Mat_C : constant Real_Matrix := Dense (Sp (3));
   Mat_D : constant Real_Matrix := Dense (Sp (4));

end Dense_AD.Integrator;
