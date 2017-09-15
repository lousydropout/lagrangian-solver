with Ada.Text_IO;

package body Auto_Differentiation is
   function "+" (X : in Real;
		 Y : in AD_Type) return AD_Type is
      N : Nat := Length (Y);
   begin
      return Const (X, N) + Y;
   end "+";
   
   function "-" (X : in Real;
		 Y : in AD_Type) return AD_Type is
      N : Nat := Length (Y);
   begin
      return Const (X, N) - Y;
   end "-";
   function "/" (X : in Real;
		 Y : in AD_Type) return AD_Type is
      N : Nat := Length (Y);
   begin
      return Const (X, N) / Y;
   end "/";

   function Const (X : in Real;
		   N : in Nat) return AD_Type is
      Y : Real_Vector (1 .. N) := (others => 0.0);
   begin
      case Level is
	 when Value =>
	    return (N, X, G0, H0);
	 when Gradient =>
	    return (N, X, Sparse (Y), H0);
	 when Hessian =>
	    return (N, X, Sparse (Y), Zero (N));
      end case;
   end Const;
   
   function Var (X    : in Real;
   		 I, N : in Nat;
   		 Dx   : in Real	:= 1.0) return AD_Type is
      Y : Real_Vector (1 .. N) := (others => 0.0);
   begin
      Y (I) := Dx;
      case Level is
	 when Value =>
	    return (N => N,
		    Val => X,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (N => N,
		    Val => X,
		    Grad => Sparse (Y),
		    Hessian => H0);
	 when Hessian =>
	    return (N => N,
		    Val => X,
		    Grad => Sparse (Y),
		    Hessian => zero (N));
      end case;
   end Var;
   
   
   function Var (X	: in Real_Vector;
   		 Length	: in Nat;
   		 Start	: in Nat := 1) return AD_Vector is
      Result : AD_Vector (1 .. X'Length);
   begin
      for I in X'Range loop
   	 Result (I - X'First + 1) := Var (X  => X (I), 
   					  I  => I - X'First + Start,
   					  N  => Length);
      end loop;
      return Result;
   end Var;
   
   function Zero (N : in Nat) return AD_Type is
      Result : AD_Type;
   begin
      Result.N := N;
      Result.Val := 0.0;
      Set_Length (Result.Grad, N);
      Result.Hessian := Zero (N);
      return Result;
   end Zero;

   
   function Val (X : in AD_Type) return Real is (X.Val);
   function Grad (X : in AD_Type) return Sparse_Vector is (X.Grad);
   --  function Grad (X : in AD_Type) return Real_Vector is (To_Array (X.Grad));
   function Hessian (X : in AD_Type) return Sparse_Matrix is (X.Hessian);
   function Length (X : in AD_Type) return Pos is (X.N);
   
   function "+" (X, Y : in AD_Type) return AD_Type is
   begin
      case Level is
	 when Value =>
	    return (N       => X.N, 
		    Val     => X.Val  + Y.Val, 
		    Grad    => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (N       => X.N, 
		    Val     => X.Val  + Y.Val, 
		    Grad    => X.Grad + Y.Grad,
		    Hessian => H0);
	 when Hessian =>
	    return (N       => X.N, 
		    Val     => X.Val  + Y.Val, 
		    Grad    => X.Grad + Y.Grad,
		    Hessian => X.Hessian + Y.Hessian);
      end case;
   end "+";
   
   function "-" (X, Y : in AD_Type) return AD_Type is
   begin
      case Level is
	 when Value =>
	    return (N       => X.N, 
		    Val     => X.Val  - Y.Val, 
		    Grad    => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (N       => X.N, 
		    Val     => X.Val  - Y.Val, 
		    Grad    => X.Grad - Y.Grad,
		    Hessian => H0);
	 when Hessian =>
	    return (N       => X.N, 
		    Val     => X.Val  - Y.Val, 
		    Grad    => X.Grad - Y.Grad,
		    Hessian => X.Hessian - Y.Hessian);
      end case;
   end "-";
   
   function "*" (X, Y : in AD_Type) return AD_Type is
      Tmp : Sparse_Matrix;
   begin
      case Level is
	 when Value =>
	    return (N       => X.N, 
		    Val     => X.Val  * Y.Val, 
		    Grad    => G0,
		    Hessian => H0);
	 when Gradient =>
	    Tmp := X.Grad * Y.Grad;
	    return (N       => X.N, 
		    Val     => X.Val  * Y.Val, 
		    Grad    => X.Val * Y.Grad + Y.Val * X.Grad,
		    Hessian => H0);
	 when Hessian =>
	    Tmp := X.Grad * Y.Grad;
	    return (N       => X.N, 
		    Val     => X.Val  * Y.Val, 
		    Grad    => X.Val * Y.Grad + Y.Val * X.Grad,
		    Hessian => 
		      X.Val * Y.Hessian + 
		      Tmp + Transpose (Tmp) +
		      Y.Val * X.Hessian);
      end case;
   end "*";
   
   function "/" (X, Y : in AD_Type) return AD_Type is
      Z  : constant Real := 1.0 / Y.Val;
      Z2 : constant Real := Z * Z;
      Z3 : constant Real := Z2 * Z;
      Tmp : Sparse_Matrix;
   begin
      case Level is
	 when Value => 
	    return (N    => X.N,
		    Val  => X.Val * Z,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient => 
	    return (N    => X.N,
		    Val  => X.Val * Z,
		    Grad => Z * X.Grad - (Z2 * X.Val) * Y.Grad,
		    Hessian => H0);
	 when Hessian => 
	    Tmp := X.Grad * Y.Grad;
	    return (N    => X.N,
		    Val  => X.Val * Z,
		    Grad => Z * X.Grad - (Z2 * X.Val) * Y.Grad,
		    Hessian => 
		      -Z2 * Tmp
		      + 2.0 * X.Val * Z3 * Y.Grad * Y.Grad 
		      - X.Val * Z2 * Y.Hessian
		      - Z2 * Transpose (Tmp)
		      + Z * X.Hessian);
      end case;
   end "/";
   
   function "**" (X : in AD_Type; N : in Integer) return AD_Type is
      Y : Real;
      Z : Real;
      H : Sparse_Matrix;
   begin
      case N is
      	 when 0 =>
      	    return (0.0 * Zero (N));
      	 when 1 =>
      	    return X;
      	 when others =>
	    case Level is
	       when Value =>
		  return (N    => X.N,
			  Val  => X.Val ** Natural (N),
			  Grad => G0,
			  Hessian => H0);
	       when Gradient =>
		  Y := Real (N) * X.Val ** Natural (N - 1);
		  return (N    => X.N,
			  Val  => X.Val ** Natural (N),
			  Grad => Y * X.Grad,
			  Hessian => H0);
	       when Hessian =>
		  Y := Real (N) * X.Val ** Natural (N - 1);
		  Z := Real (N * (N - 1)) * X.Val ** Natural (N - 2);
		  H := Y * X.Hessian + Z * X.grad * X.grad;
		  return (N    => X.N,
			  Val  => X.Val ** Natural (N),
			  Grad => Y * X.Grad,
			  Hessian => H);
	    end case;
      end case;
   end "**";

   
   
   function Sin (X : in AD_Type) return AD_Type is
      use Real_Functions;
      S : constant Real := Sin (X.Val);
      C : constant Real := Cos (X.Val);
   begin
      case Level is
	 when Value =>
	    return (N => X.N,
		    Val => S,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (N => X.N,
		    Val => S,
		    Grad => C * X.Grad,
		    Hessian => H0);
	 when Hessian =>
	    return (N => X.N,
		    Val => S,
		    Grad => C * X.Grad,
		    Hessian => -S * X.Grad * X.Grad + C * X.Hessian);
      end case;
   end Sin;
   function Cos (X : in AD_Type) return AD_Type is
      use Real_Functions;
      S : constant Real := Sin (X.Val);
      C : constant Real := Cos (X.Val);
   begin
      case Level is
	 when Value =>
	    return (N    => X.N, 
		    Val  => C,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (N    =>  X.N, 
		    Val  =>  C,
		    Grad => -S * X.Grad,
		    Hessian => H0);
	 when Hessian =>
	    return (N    =>  X.N, 
		    Val  =>  C,
		    Grad => -S * X.Grad,
		    Hessian => -C * X.Grad * X.Grad - S * X.Hessian);
      end case;
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
      case Level is
	 when Value =>
	    return (N    => X.N,
		    Val  => T,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    Z := Y * X.Grad;
	    return (N    => X.N,
		    Val  => T,
		    Grad => Z,
		    Hessian => H0);
	 when Hessian =>
	    Z := Y * X.Grad;
	    return (N    => X.N,
		    Val  => T,
		    Grad => Z,
		    Hessian => (1.0 + Y) * X.Hessian + (2.0 * T) * Z * X.Grad);
      end case;
   end Tan;
	      
      
   function Exp (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : constant Real := Exp (X.Val);
      G : Sparse_Vector;
   begin
      case Level is
	 when Value =>
	    return (N    => X.N,
		    Val  => Y,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    G := Y * X.Grad;
	    return (N    => X.N,
		    Val  => Y,
		    Grad => G,
		    Hessian => H0);
	 when Hessian =>
	    G := Y * X.Grad;
	    return (N    => X.N,
		    Val  => Y,
		    Grad => G,
		    Hessian => Y * X.Hessian + G * X.Grad);
      end case;
   end Exp;
   
   
   function Log (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : Real;
      Z : Sparse_Vector;
   begin
      pragma Assert (X.Val > 0.0);
      Y := 1.0 / X.Val;
      case Level is
	 when Value =>
	    return (N    => X.N,
		    Val  => Log (X.Val),
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    Z := Y * X.Grad;
	    return (N    => X.N,
		    Val  => Log (X.Val),
		    Grad => Z,
		    Hessian => H0);
	 when Hessian =>
	    Z := Y * X.Grad;
	    return (N    => X.N,
		    Val  => Log (X.Val),
		    Grad => Z,
		    Hessian => -Z * Z + Y * X.Hessian);
      end case; 
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
   
   
   
   function "+" (X, Y : in AD_2D) return AD_2D is
   begin
      return (X (1) + Y (1), 
	      X (2) + Y (2));
   end "+";
   
   function "-" (X, Y : in AD_2D) return AD_2D is
   begin
      return (X (1) - Y (1), 
	      X (2) - Y (2));
   end "-";
   
   function Dot (X, Y : in AD_2D) return AD_Type is
   begin
      return (X (1) * Y (1) + X (2) * Y (2));
   end Dot;



   procedure Print (X : in AD_Type) is
      use Ada.Text_IO;
   begin
      Put_Line ("   value: " & Real'Image (X.Val));
      Print (X.Grad);
      New_Line;
   end Print;
   
end Auto_Differentiation;
