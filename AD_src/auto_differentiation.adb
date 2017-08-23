with Ada.Text_IO;
package body Auto_Differentiation is
   
   function Plus1 (X : in Real) return Real is
   begin
      return (X + 1.0);
   end Plus1;
   
   
   function Var (X    : in Real;
   		 I, N : in Nat;
   		 Dx   : in Real	:= 1.0) return AD_Type is
      Y : Real_Array (1 .. N) := (others => 0.0);
   begin
      Y (I) := Dx;
      return (N => N,
   	      Val => X,
   	      Grad => Sparse (Y),
   	      Hessian => zero (N));
   end Var;
   
   
   --  function Var (X	: in Real_Array;
   --  		 Length	: in Nat;
   --  		 Start	: in Nat := 1) return AD_Vector is
   --     Result : AD_Vector (1 .. X'Length);
   --  begin
   --     for I in X'Range loop
   --  	 Result (I - X'First + 1) := Var (X  => X (I), 
   --  					  I  => I - X'First + Start,
   --  					  N  => Length);
   --     end loop;
   --     return Result;
   --  end Var;
   
   function Val (X : in AD_Type) return Real is (X.Val);
   function Grad (X : in AD_Type) return Sparse_Vector is (X.Grad);
   function Grad (X : in AD_Type) return Real_Array is (To_Array (X.Grad));
   function Hessian (X : in AD_Type) return Sparse_Matrix is (X.Hessian);
   function Length (X : in AD_Type) return Pos is (X.N);
   
   function "+" (X, Y : in AD_Type) return AD_Type is
   begin
      return (N       => X.N, 
   	      Val     => X.Val  + Y.Val, 
   	      Grad    => X.Grad + Y.Grad,
   	      Hessian => X.Hessian + Y.Hessian);
   end "+";
   
   function "-" (X, Y : in AD_Type) return AD_Type is
   begin
      return (N       => X.N, 
   	      Val     => X.Val  - Y.Val, 
   	      Grad    => X.Grad - Y.Grad,
   	      Hessian => X.Hessian - Y.Hessian);
   end "-";
   
   function "*" (X, Y : in AD_Type) return AD_Type is
      Tmp : Sparse_Matrix := X.Grad * Y.Grad;
   begin
      return (N       => X.N, 
   	      Val     => X.Val  * Y.Val, 
   	      Grad    => X.Val * Y.Grad + Y.Val * X.Grad,
   	      Hessian => 
   		X.Val * Y.Hessian +
   		Tmp + Transpose (Tmp) +
   		Y.Val * X.Hessian);
   end "*";
   
   function "/" (X, Y : in AD_Type) return AD_Type is
      Z : constant Real := 1.0 / Y.Val;
      Tmp : constant Sparse_Matrix := X.Grad * Y.Grad;
      Z2 : constant Real := Z * Z;
      Z3 : constant Real := Z2 * Z;
   begin
      return (N    => X.N,
   	      Val  => X.Val * Z,
   	      Grad => Z * X.Grad - (Z ** 2 * X.Val) * Y.Grad,
   	      Hessian => 
   		-Z2 * Tmp
   		+ 2.0 * X.Val * Z3 * Y.Grad * Y.Grad 
   		- X.Val * Z2 * Y.Hessian
   		- Z2 * Transpose (Tmp)
   		+ Z * X.Hessian);
   end "/";
   
   function "**" (X : in AD_Type; N : Pos) return AD_Type is
      Y : constant Real := Real (N) * X.Val ** Natural (N - 1);
      Z : constant Real := Real (N * (N - 1)) * X.Val ** Natural (N - 2);
   begin
      return (N    => X.N,
   	      Val  => X.Val ** Natural (N),
   	      Grad => Y * X.Grad,
   	      Hessian => Y * X.Hessian + Z * X.Grad * X.Grad);
   end "**";

   
   
   function Sin (X : in AD_Type) return AD_Type is
      use Real_Functions;
      S : constant Real := Sin (X.Val);
      C : constant Real := Cos (X.Val);
   begin
      return (N => X.N,
   	      Val => S,
   	      Grad => C * X.Grad,
   	      Hessian => -S * X.Grad * X.Grad + C * X.Hessian);
   end Sin;
   function Cos (X : in AD_Type) return AD_Type is
      use Real_Functions;
      S : constant Real := Sin (X.Val);
      C : constant Real := Cos (X.Val);
   begin
      return (N    =>  X.N, 
   	      Val  =>  C,
   	      Grad => -S * X.Grad,
   	      Hessian => -C * X.Grad * X.Grad - S * X.Hessian);
   end Cos;
   function Tan (X : in AD_Type) return AD_Type is
      use Real_Functions;
      C : constant Real := Cos (X.Val);
      T : Real;
      Y : Real;
      Z : Sparse_Vector;
   begin
      pragma Assert (C /= 0.0);
      T := Sin (X.Val) / C;
      Y := (1.0 - T ** 2);
      Z := Y * X.Grad;
      return (N    => X.N,
   	      Val  => T,
   	      Grad => Z,
   	      Hessian => (1.0 + Y) * X.Hessian + (2.0 * T) * Z * X.Grad);
   end Tan;
	      
      
   function Exp (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : constant Real := Exp (X.Val);
      G : constant Sparse_Vector := Y * X.Grad;
   begin
      return (N    => X.N,
   	      Val  => Y,
   	      Grad => G,
   	      Hessian => Y * X.Hessian + G * X.Grad);
   end Exp;
   
   
   function Log (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : Real;
      Z : Sparse_Vector;
   begin
      pragma Assert (X.Val > 0.0);
      Y := 1.0 / X.Val;
      Z := Y * X.Grad;
      return (N    => X.N,
   	      Val  => Log (X.Val),
   	      Grad => Z,
   	      Hessian => -Z * Z + Y * X.Hessian);
   end Log;
   
   function "-" (X : in AD_Type) return AD_Type is
   begin
      return (N    =>  X.N, 
   	      Val  => -X.Val,
   	      Grad => -X.Grad,
   	      Hessian => -X.Hessian);
   end "-";
   
      
   function "*" (Y : in Real; X : in AD_Type) return AD_Type is
   begin
      return (N    => X.N,
   	      Val  => Y * X.Val,
   	      Grad => Y * X.Grad,
   	      Hessian => Y * X.Hessian);
   end "*";
   


   procedure Print (X : in AD_Type) is
      use Ada.Text_IO;
   begin
      Put_Line ("   value: " & Real'Image (X.Val));
      Print (X.Grad);
      New_Line;
   end Print;
   
end Auto_Differentiation;