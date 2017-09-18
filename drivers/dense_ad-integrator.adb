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
      
      Mat_C := D and Eye2N; -- temporary, overwritten below
      Mat_A := Mat_C * (EyeK and TL);
      Mat_B := Mat_C * (EyeK and BR);
      Mat_C := EyeK and TR;
      Mat_D := EyeK and BL;
   end Setup;
   
   
   function Collocation (Lagrangian : not null access 
			   function (X : Vector) return AD_Type;
			 Var        : in     Variable;
			 Control    : in out Control_Type) return Real_Vector is
      use Ada.Text_IO, Real_IO;
      Old : constant Evaluation_Level :=  Level;
      
      Q   : Real_Vector (1 .. NK);
      DQ  : Sparse_Vector;
      F   : Sparse_Vector;
      J   : Sparse_Matrix;
      Res : Real := 1.0;
   begin
      Level := Hessian;
      ------------------------------------------------
      for I in 1 .. K loop
	 Q ((I - 1) * N + 1 .. I * N) := Var.X;
      end loop;
      ------------------------------------------------
      while Res > 1.0e-10 loop
	 FJ (Lagrangian, Var, Control, Q, F, J);
	 DQ := Numerics.Sparse_Matrices.CSparse.Solve (J, F);
	 Q (N + 1 .. N * K) := Q (N + 1 .. N * K) - To_Array (DQ);
	 Res := Norm (F);
	 --  Put (Res); New_Line;
      end loop;
      --  New_Line;
      ------------------------------------------------
      Control.Err := Norm (DQ) / Norm (Q);
      Level := Old;

      return Q;
   end Collocation;

   
   procedure FJ (Lagrangian : not null access 
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
      A    : constant Sparse_Matrix := Mat_A - Control.Dt * Mat_C;
      B    : constant Sparse_Matrix := Mat_B - Control.Dt * Mat_D;
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
      F := A * Q - B * U; F := Remove_1stN (F, N);
      J := A     - B * V; J := Remove_1stN (J, N);
   end FJ;
   

begin
   
   Setup;
   
end Dense_AD.Integrator;
