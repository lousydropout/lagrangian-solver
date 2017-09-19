with Ada.Text_IO, Numerics.Dense_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Ada.Text_IO, Numerics.Dense_Matrices;
package body Dense_AD.Integrator is
   
   procedure Print_Lagrangian (Var	  : in     Variable;
			       Lagrangian : not null access
				 function (X : Vector) return AD_Type;
			       Fore : in Field := 3;
			       Aft  : in Field := 5;
			       Exp  : in Field := 3) is
   begin
      Print_Lagrangian (File => Standard_Output, 
			Var => Var, 
			Lagrangian => Lagrangian,
			Fore => Fore,
			Aft => Aft,
			Exp => Exp);
   end Print_Lagrangian;

   
   procedure Print_Lagrangian (File	  : in     Ada.Text_IO.File_Type;
			       Var	  : in     Variable;
			       Lagrangian : not null access
				 function (X : Vector) return AD_Type;
			       Fore : in Field := 3;
			       Aft  : in Field := 5;
			       Exp  : in Field := 3) is
      use Real_IO;
      Old : constant Evaluation_Level := Level;
      L   : AD_Type;
      P   : Real_Vector (Half_N + 1 .. N);
      V   : Real_Vector renames Var.X (Half_N + 1 .. N);
      H   : Real;
   begin
      Level := Gradient;
      L := Lagrangian (Var.X);
      P := Grad (L) (Half_N + 1 .. N);
      H := Dot (V, P) - Val (L);
      
      -- Print time
      Put (File, Var.T, Fore, Aft, Exp); Put (File, ",  ");
      -- Print generalized positions and velocities
      for I in 1 .. N loop
	 Put (File, Var.X (I), Fore, Aft, Exp); Put (File, ",  ");
      end loop;
      -- Print canonical momenta
      for I in Half_N + 1 .. N loop
	 Put (File, P (I), Fore, Aft, Exp); Put (File, ",  ");
      end loop;
      -- Print Hamiltonian
      Put (File, H, Fore => Fore, Exp => Exp); New_Line (File);
      
      Level := Old;
   end Print_Lagrangian;
   
   
   function Interpolate (A : in Array_Of_Vectors;
			 T : in Real;
			 L : in Real;
			 R : in Real) return Vector is
      X : Vector;
   begin
      for I in 1 .. N loop
	 X (I) := Interpolate (A (I), T, L, R);
      end loop;
      return X;
   end Interpolate;
   
   
   procedure Update (Var : in out Variable;
		     Y	 : in     Real_Vector;
		     Dt	 : in     Real) is
   begin
      Var.X := Y (Y'Last - N + 1 .. Y'Last);
      Var.T := Var.T + Dt;
   end Update;
   
   function Update (Lagrangian : not null access 
		      function (X : Vector) return AD_Type;
		    Var        : in     Variable;
		    Control    : in out Control_Type;
		    Density    : in     Dense_Or_Sparse) return Real_Vector is
      use Real_Functions;
      Y      : Real_Vector (1 .. NK);
      A      : Array_Of_Vectors;
      Dt     : Real;
      Time   : Real;
   begin
      Collocation_Density := Density;
      Control.Err := 1.0;
      Control.Dt  := Control.Dtn; 
      -- Set initial guess for Y -----------------------------------------
      for I in 1 .. K loop Y ((I - 1) * N + 1 .. I * N) := Var.X; end loop;
      ---------------------------------------------------------------------
      while Control.Err > Control.Eps loop
	 Iterate (Lagrangian, Y, Var, Control);
	 Dt := 0.8 * Control.Dt
	   * (Control.Eps / (Control.Err + 1.0e-40)) ** (1.0 / Real (K - 1));
	 if Control.Err > Control.Eps then
	    -- Since the error is above the threshold, redo the calculation
	    --   for a smaller Dt. But use above estimate of Y for initial
	    --   guess.
	    A := Chebyshev_Transform (Y);
	    for I in 1 .. K loop
	       Time := Var.T + Dt * Grid (I);
	       Y ((I - 1) * N + 1 .. I * N)
		 := Interpolate (A, Time, Var.T, Var.T + Control.Dt);
	    end loop;
	    Control.Dt := Dt;
	 end if;
      end loop;
      Control.Dtn := Dt; -- Set future Dt
      return Y;
   end Update;
   
   function Split (Y : in     Real_Vector) return Array_Of_Vectors is
      A : Array_Of_Vectors;
   begin
      for J in 1 .. N loop
      	 for I in 1 .. K loop
      	    A (J) (I) := Y (N * (I - 1) + J);
	 end loop;
      end loop;
      return A;
   end Split;
   
   
   function Chebyshev_Transform (Y : in Real_Vector) return Array_Of_Vectors is
      A : Array_Of_Vectors := Split (Y);
   begin
      for J in 1 .. N loop
	 A (J) := CGL_Transform (A (J));
      end loop;
      return A;
   end Chebyshev_Transform;
      
   procedure Colloc (Lagrangian : not null access 
		       function (X : Vector) return AD_Type;
		     Q          : in out Real_Vector;
		     Var        : in     Variable;
		     Control    : in out Control_Type) is
   begin
      case Collocation_Density is
	 when Sparse => Sp_Collocation (Lagrangian, Q, Var, Control);
	 when Dense =>     Collocation (Lagrangian, Q, Var, Control);
      end case;
   end Colloc;
   
   procedure Iterate (Lagrangian : not null access 
			function (X : Vector) return AD_Type;
		      Y          : in out Real_Vector;
		      Var        : in     Variable;
		      Control    : in out Control_Type) is
      Y1, Y2 : Real_Vector (1 .. NK);
      Var2 : Variable;
      C2   : Control_Type := Control;
      A1   : Array_Of_Vectors;
      A2   : Array_Of_Vectors;
      Err  : Real := 1.0;
      Time : Real;
      Tf   : constant Real := Var.T + Control.Dt;
   begin
      -------------------------------------------------------
      C2.Dt := 0.5 * C2.Dt;
      -------------------------------------------------------
      -- Interpolate to get initial guess for Y1
      A1 := Chebyshev_Transform (Y);
      for I in 1 .. K loop
      	 Time := Var.T + C2.Dt * Grid (I);
      	 Y1 ((I - 1) * N + 1 .. I * N)
      	   := Interpolate (A1, Time, Var.T, Var.T + Control.Dt);
      end loop;
      -------------------------------------------------------
      -- Iterate from Y1 to Y2
      Colloc (Lagrangian, Y1, Var, C2);
      Var2.X := Y1 (NK - N + 1 .. NK);
      Var2.T := Var.T + C2.Dt;
      -------------------------------------------------------
      -- Set intial guess for Y2
      for I in 1 .. K loop Y2 ((I - 1) * N + 1 .. I * N) := Var2.X; end loop;
      -------------------------------------------------------
      -- Update Y2
      Colloc (Lagrangian, Y2, Var2, C2);
      -------------------------------------------------------
      -- Interpolate to get initial guess for Y
      A1 := Chebyshev_Transform (Y1);
      A2 := Chebyshev_Transform (Y2);
      -------------------------------------------------------
      for I in 1 .. K loop
	 Time := Var.T + Control.Dt * Grid (I);
	 if Time <= Var.T + C2.Dt then
	    Y1 ((I - 1) * N + 1 .. I * N) 
	      := Interpolate (A1, Time, Var.T, Var2.T);
	 else
	    Y1 ((I - 1) * N + 1 .. I * N)
	      := Interpolate (A2, Time, Var2.T, Tf);
	 end if;
      end loop;
      -------------------------------------------------------
      -- Update Y
      Y := Y1;
      Colloc (Lagrangian, Y1, Var, Control);
      -------------------------------------------------------
      Control.Err := Norm (Y - Y1);
   end Iterate;
   
   function Setup return Array_Of_Sparse_Matrix is
      Sp  : Array_Of_Sparse_Matrix;
      Tmp : Sparse_Matrix;
   begin
      Tmp := D and Eye2N; -- temporary, overwritten below
      Sp (1) := Tmp * (EyeK and Top_Left);
      Sp (2) := Tmp * (EyeK and Bottom_Right);
      Sp (3) := EyeK and Top_Right;
      Sp (4) := EyeK and Bottom_Left;
      return Sp;
   end Setup;
   
   
   procedure Collocation (Lagrangian : not null access 
			   function (X : Vector) return AD_Type;
			 Q          : in out Real_Vector;
			 Var        : in     Variable;
			 Control    : in out Control_Type) is
      --  use Ada.Text_IO, Real_IO;
      Old : constant Evaluation_Level :=  Level;
      
      DQ  : Real_Vector (1 .. NK - N);
      F   : Real_Vector (1 .. NK - N);
      J   : Real_Matrix (1 .. NK - N, 1 .. NK - N);
      Res : Real := 1.0;
      It  : Pos  := 0;
   begin
      Level := Hessian;
      ------------------------------------------------
      while Res > 1.0e-10 loop
	 It := It + 1;
	 FJ (Lagrangian, Var, Control, Q, F, J);
	 DQ := Solve (J, F);
	 Q (N + 1 .. N * K) := Q (N + 1 .. N * K) - DQ;
	 Res := Norm (F);
	 --  Put (Res); New_Line;
	 if It > 10 then exit; end if;
      end loop;
      --  New_Line;
      ------------------------------------------------
      Level := Old;
   end Collocation;

   
   procedure FJ (Lagrangian : not null access 
		    function (X : Vector) return AD_Type;
		  Var     : in     Variable;
		  Control : in     Control_Type;
		  Q       : in     Real_Vector;
		  F       :    out Real_Vector;
		  J       :    out Real_Matrix) is
      L    : AD_Type;
      X    : Vector;
      Tmp  : Integer;
      U, R : Real_Vector (1 .. NK) := (others => 0.0);
      V, S : Real_Matrix (1 .. NK, 1 .. NK) := (others => (others => 0.0));
      A    : constant Real_Matrix := Mat_A - Control.Dt * Mat_C;
      B    : constant Real_Matrix := Mat_B - Control.Dt * Mat_D;
   begin
      pragma Assert (Q'Length = N * K);
      -------------------------------------
      for I in 1 .. K loop
	 Tmp := N * (I - 1);
	 X := Q (Tmp + 1 .. Tmp + N);
	 L := Lagrangian (X);
	 Copy (From => Grad (L), To => U, Start => Tmp + 1);
	 Copy (From => Hessian (L), 
	       To => V, 
	       Start_I => Tmp + 1, 
	       Start_J => Tmp + 1);
      end loop;
      -------------------------------------
      R := A * Q - B * U; F := Remove_1st_N (R, N);
      S := A     - B * V; J := Remove_1st_N (S, N);
   end FJ;
   
   
   procedure Sp_Collocation (Lagrangian : not null access 
			       function (X : Vector) return AD_Type;
			     Q          : in out Real_Vector;
			     Var        : in     Variable;
			     Control    : in out Control_Type) is
      --  use Ada.Text_IO, Real_IO;
      Old : constant Evaluation_Level :=  Level;
      
      DQ  : Sparse_Vector;
      F   : Sparse_Vector;
      J   : Sparse_Matrix;
      Res : Real := 1.0;
      It  : Pos  := 0;
   begin
      Level := Hessian;
      ------------------------------------------------
      while Res > 1.0e-10 loop
	 It := It + 1;
   	 Sp_FJ (Lagrangian, Var, Control, Q, F, J);
   	 DQ := Numerics.Sparse_Matrices.CSparse.Solve (J, F);
   	 Q (N + 1 .. N * K) := Q (N + 1 .. N * K) - To_Array (DQ);
   	 Res := Norm (F);
   	 --  Put (Res); New_Line;
	 if It > 10 then exit; end if;
      end loop;
      --  New_Line;
      ------------------------------------------------
      Level := Old;

   end Sp_Collocation;

   
   procedure Sp_FJ (Lagrangian : not null access 
		      function (X : Vector) return AD_Type;
		    Var     : in     Variable;
		    Control : in     Control_Type;
		    Q       : in     Real_Vector;
		    F       :    out Sparse_Vector;
		    J       :    out Sparse_Matrix) is
      L    : AD_Type;
      X    : Real_Vector (1 .. N);
      Tmp  : Integer;
      U    : Sparse_Vector;
      V    : Sparse_Matrix;
      A    : constant Sparse_Matrix := Sp (1) - Control.Dt * Sp (3);
      B    : constant Sparse_Matrix := Sp (2) - Control.Dt * Sp (4);
   begin
      pragma Assert (Q'Length = N * K);
      -------------------------------------
      X := Q (1 .. N);
      L := Lagrangian (X);
      U := Sparse (Grad (L));
      V := Sparse (Hessian (L));
      for I in 2 .. K loop
   	 Tmp := N * (I - 1);
   	 X := Q (Tmp + 1 .. Tmp + N);
   	 L := Lagrangian (X);
   	 U := U or Sparse (Grad (L));
   	 V := V or Sparse (Hessian (L));
      end loop;
      -------------------------------------
      F := A * Q - B * U; F := Remove_1st_N (F, N);
      J := A     - B * V; J := Remove_1st_N (J, N);
   end Sp_FJ;
   
   
   
   procedure Print_XYZ (File : in File_Type;
			Var  : in Variable) is
      use Real_Functions, Real_IO;
      X :  Real_Vector renames Var.X;
      X1, Y1, X2, Y2 : Real;
      R : constant Real := 10.0;
   begin
      X1 := -R * Sin (X (1));
      Y1 :=  R * Cos (X (1));
      X2 := X1 - R * Sin (2.0 * X (1) + X (2));
      Y2 := Y1 + R * Cos (2.0 * X (1) + X (2));
      -- print header
      Put_Line (File, "3");
      Put (File, "Properties=pos:R:2   Time=");
      Put (File, Var.T, Aft => 5, Exp => 0);
      New_Line (File);
      -- position of ball 1
      Put_Line (File, "0.0     0.0     5.0");
      -- position of ball 2
      Put (File => File, Item => X1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y1);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
      -- position of ball 3
      Put (File => File, Item => X2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => Y2);
      Put (File => File, Item => "     ");
      Put (File => File, Item => "5.0");
      New_Line (File => File);
   end Print_XYZ;

   
end Dense_AD.Integrator;
