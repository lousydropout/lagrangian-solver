with Ada.Text_IO, Numerics.Dense_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Ada.Text_IO, Numerics.Dense_Matrices;
package body Dense_AD.Integrator is
   
   procedure Print_Lagrangian (X : in Vector;
			       Lagrangian : not null access
				 function (X : Vector)  return AD_Type) is
      use Real_IO;
      L : AD_Type := Lagrangian (X);
   begin
      Put ("Val (L) = ");
      Put (Val (L), Aft => 3, Exp => 0);
      New_Line;
      Put_Line ("Grad (L) = ");
      Print (Grad (L));
      Put_Line ("Hessian (L) = ");
      Print (Hessian (L));
   end Print_Lagrangian;
   
   
   function Update (Lagrangian : not null access 
		      function (X : Vector) return AD_Type;
		    Var        : in     Variable;
		    Control    : in out Control_Type;
		    Choose     : in     Dense_Or_Sparse) return Real_Vector is
      use Real_Functions;
      Y    : Real_Vector (1 .. NK);
      A    : array (1 .. 2 * N) of Real_Vector (1 .. K);
      Dt   : Real;
      Time : Real;
   begin
      Control.Err := 1.0;
      Control.Dt  := Control.Dtn;
      -- Set initial guess for Y
      for I in 1 .. K loop Y ((I - 1) * N + 1 .. I * N) := Var.X; end loop;
      while Control.Err > Control.Eps loop
	 Iterate (Lagrangian, Y, Var, Control, Choose);
	 Dt := 0.8 * Control.Dt * (Control.Eps / (Control.Err + 1.0e-40)) ** (1.0 / Real (K - 1));
	 if Control.Err > Control.Eps then
	    for J in 1 .. N loop
	       for I in 1 .. K loop
		  A (J) (I) := Y (N * (I - 1) + J);
	       end loop;
	    end loop;
	    for J in 1 .. N loop
	       A (J) := CGL_Transform (A (J));
	    end loop;
	    for I in 1 .. K loop
	       Time := Var.T + Dt * Grid (I);
	       for J in 1 .. N loop
		  Y ((I - 1) * N + J)
		    := Interpolate (A (J), Time, Var.T, Var.T + Control.Dt);
	       end loop;
	    end loop;
	    Control.Dt := Dt;
	 end if;
      end loop;
      Control.Dtn := Dt;
      return Y;
   end Update;
   
   
   
   procedure Iterate (Lagrangian : not null access 
			function (X : Vector) return AD_Type;
		      Y          : in out Real_Vector;
		      Var        : in     Variable;
		      Control    : in out Control_Type;
		      Choose     : in     Dense_Or_Sparse) is
      Y1, Y2 : Real_Vector (1 .. NK);
      Var2 : Variable;
      C2   : Control_Type := Control;
      A1   : array (1 .. 2 * N) of Real_Vector (1 .. K);
      A2   : array (1 .. 2 * N) of Real_Vector (1 .. K);
      Err  : Real := 1.0;
      Time : Real;
      DT   : constant Real := Var.T + Control.Dt;
   begin
      C2.Dt := 0.5 * C2.Dt;
      -------------------------------------------------------
      -- Interpolate to get initial guess for Y1
      for J in 1 .. N loop
      	 for I in 1 .. K loop
      	    A1 (J) (I) := Y (N * (I - 1) + J);
	 end loop;
      end loop;
      for J in 1 .. N loop
      	 A1 (J) := CGL_Transform (A1 (J));
      end loop;
      for I in 1 .. K loop
	 Time := Var.T + C2.Dt * Grid (I);
	 for J in 1 .. N loop
	    Y1 ((I - 1) * N + J)
	      := Interpolate (A1 (J), Time, Var.T, Var.T + C2.Dt);
	 end loop;
      end loop;
      -------------------------------------------------------
      -- Iterate from Y1 to Y2
      for I in 1 .. K loop Y1 ((I - 1) * N + 1 .. I * N) := Var.X; end loop;
      case Choose is
	 when Sparse => Sp_Collocation (Lagrangian, Y1, Var, C2);
	 when Dense =>     Collocation (Lagrangian, Y1, Var, C2);
      end case;
      Var2.X := Y1 (NK - N + 1 .. NK);
      Var2.T := Var.T + C2.Dt;
      -------------------------------------------------------
      -- Interpolate to get intial guess for Y2
      for I in 1 .. K loop 
	 Time := Var2.T + C2.Dt * Grid (I);
	 for J in 1 .. N loop
	    Y2 ((I - 1) * N + J)
	      := (1.0 - Grid (I)) * Var2.X (J) 
	      + Grid (I) * Interpolate (A1 (J), Time, Var2.T, DT);
	 end loop;
      end loop;
      -------------------------------------------------------
      -- Update Y2
      case Choose is
	 when Sparse => Sp_Collocation (Lagrangian, Y2, Var2, C2);
	 when Dense =>     Collocation (Lagrangian, Y2, Var2, C2);
      end case;
      -------------------------------------------------------
      -- Interpolate to get initial guess for Y
      for J in 1 .. N loop
      	 for I in 1 .. K loop
      	    A1 (J) (I) := Y1 (N * (I - 1) + J);
      	    A2 (J) (I) := Y2 (N * (I - 1) + J);
      	 end loop;
      end loop;
      -------------------------------------------------------
      for J in 1 .. N loop
      	 A1 (J) := CGL_Transform (A1 (J));
      	 A2 (J) := CGL_Transform (A2 (J));
      end loop;
      -------------------------------------------------------
      for I in 1 .. K loop
	 Time := Var.T + Control.Dt * Grid (I);
	 if Time <= Var.T + C2.Dt then
	    for J in 1 .. N loop
	       Y1 ((I - 1) * N + J) := Interpolate (A1 (J), Time, Var.T, Var2.T);
	    end loop;
	 else
	    for J in 1 .. N loop
	       Y1 ((I - 1) * N + J) := Interpolate (A2 (J), Time, Var2.T, DT);
	    end loop;
	 end if;
      end loop;
      -------------------------------------------------------
      -- Update Y
      Y := Y1;
      case Choose is
	 when Sparse => Sp_Collocation (Lagrangian, Y1, Var, Control);
	 when Dense =>     Collocation (Lagrangian, Y1, Var, Control);
      end case;
      -------------------------------------------------------
      Control.Err := Norm (Y - Y1);
   end Iterate;
   
   procedure Setup is
      TL, TR, BL, BR : Sparse_Matrix;
      EyeN  : constant Sparse_Matrix := Eye (Half_N);
      EyeK  : constant Sparse_Matrix := Eye (K);
      Eye2N : constant Sparse_Matrix := Eye (N);
      D     : constant Sparse_Matrix := Sparse (Der);
   begin
      TL := Sparse (Top_Left)     and EyeN;
      TR := Sparse (Top_Right)    and EyeN;
      BL := Sparse (Bottom_Left)  and EyeN;
      BR := Sparse (Bottom_Right) and EyeN;
      
      Sp_C := D and Eye2N; -- temporary, overwritten below
      Sp_A := Sp_C * (EyeK and TL);
      Sp_B := Sp_C * (EyeK and BR);
      Sp_C := EyeK and TR;
      Sp_D := EyeK and BL;
      
      Dense (Sp => Sp_A, A => Mat_A);
      Dense (Sp => Sp_B, A => Mat_B);
      Dense (Sp => Sp_C, A => Mat_C);
      Dense (Sp => Sp_D, A => Mat_D);
   end Setup;
   
   
   procedure Collocation (Lagrangian : not null access 
			   function (X : Vector) return AD_Type;
			 Q          : in out Real_Vector;
			 Var        : in     Variable;
			 Control    : in out Control_Type) is
      --  use Ada.Text_IO, Real_IO;
      Old : constant Evaluation_Level :=  Level;
      
      --  Q   : Real_Vector (1 .. NK);
      DQ  : Real_Vector (1 .. NK - N);
      F   : Real_Vector (1 .. NK - N);
      J   : Real_Matrix (1 .. NK - N, 1 .. NK - N);
      Res : Real := 1.0;
      It  : Pos  := 0;
   begin
      Level := Hessian;
      ------------------------------------------------
      --  for I in 1 .. K loop
      --  	 Q ((I - 1) * N + 1 .. I * N) := Var.X;
      --  end loop;
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
      A    : constant Sparse_Matrix := Sp_A - Control.Dt * Sp_C;
      B    : constant Sparse_Matrix := Sp_B - Control.Dt * Sp_D;
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
   

begin
   
   Setup;
   
end Dense_AD.Integrator;
