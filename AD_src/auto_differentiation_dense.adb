package body Auto_Differentiation_Dense is
   
   function Var (X  : in Real;
		 I  : in Nat;
		 Dx : in Real := 1.0) return AD_Type is
      Y : Real_Vector := G0;
   begin
      Y (I) := Dx;
      case Level is
	 when Value    => return (X, G0, H0);
	 when Gradient => return (X, Y, H0);
	 when Hessian  => return (X, Y, H0);
      end case;
   end Var;
   
   function Var (X : in Real_Vector) return AD_Vector is
      Result : AD_Vector;
   begin
      for I in X'Range loop
	 Result (I - X'First + 1) := Var (X => X (I),
					  I => I - X'First + 1);
      end loop;
      return Result;
   end Var;
   
   function Const (X : in Real) return AD_Type is
   begin
      return (X, G0, H0);
   end Const;
   
   function Val (X : in AD_Type) return Real is
   begin
      return X.Val;
   end Val;
   
   function Grad (X : in AD_Type) return Real_Vector is
   begin
      return X.Grad;
   end Grad;
   
   function Hessian (X : in AD_Type) return Real_Matrix is
   begin
      return X.Hessian;
   end Hessian;
   
   function "+" (X, Y : in AD_Type) return AD_Type is
   begin
      case Level is
	 when Value => return (X.Val + Y.Val, G0, H0);
	 when Gradient => return (X.Val + Y.Val, X.Grad + Y.grad, H0);
	 when Hessian => return (X.Val + Y.Val, X.Grad + Y.grad, 
				 X.Hessian + Y.Hessian);
      end case;
   end "+";
   
   function "-" (X, Y : in AD_Type) return AD_Type is
   begin
      case Level is
	 when Value => return (X.Val - Y.Val, G0, H0);
	 when Gradient => return (X.Val - Y.Val, X.Grad - Y.grad, H0);
	 when Hessian => return (X.Val - Y.Val, X.Grad - Y.grad, 
				 X.Hessian - Y.Hessian);
      end case;
   end "-";
   
   function "*" (X, Y : in AD_Type) return AD_Type is
      Tmp : Real_Matrix (1 .. N, 1 .. N);
   begin
      case Level is
	 when Value => return (X.Val * Y.Val, G0, H0);
	 when Gradient => return (X.Val * Y.Val, 
				  Y.Val * X.Grad + X.Val * Y.grad, 
				  H0);
	 when Hessian => 
	    Tmp := Outer (X.Grad, Y.Grad);
	    return (Val => X.Val * Y.Val, 
		    Grad => X.Val * Y.Grad + Y.Val * X.Grad,
		    Hessian => X.Val * Y.Hessian + Y.Val * X.Hessian
		      + Tmp + Transpose (Tmp));
      end case;
   end "*";
     
   function "/" (X, Y : in AD_Type) return AD_Type is
      Z  : constant Real := 1.0 / Y.Val;
      Z2 : constant Real := Z * Z;
      Z3 : constant Real := Z2 * Z;
      Tmp : Real_Matrix (1 .. N, 1 .. N);
   begin
      case Level is
	 when Value => 
	    return (Val => X.Val * Z, Grad => G0, Hessian => H0);
	 when Gradient => 
	    return (Val  => X.Val * Z,
		    Grad => Z * X.Grad - (Z2 * X.Val) * Y.Grad,
		    Hessian => H0);
	 when Hessian => 
	    Tmp := Outer (X.Grad, Y.Grad);
	    return (Val  => X.Val * Z,
		    Grad => Z * X.Grad - (Z2 * X.Val) * Y.Grad,
		    Hessian => 
		      -Z2 * Tmp
		      + 2.0 * X.Val * Z3 * Y.Grad * Y.Grad 
		      - X.Val * Z2 * Y.Hessian
		      - Z2 * Transpose (Tmp)
		      + Z * X.Hessian);
      end case;
   end "/";

   function "**" (X : in AD_Type; K : in Integer) return AD_Type is
      Y : Real;
      Z : Real;
      H : Real_Matrix (1 .. N, 1 .. N);
   begin
      if K < 0 then pragma Assert (X.Val /= 0.0); end if;
      case K is
      	 when 0 => return (1.0, G0, H0);
      	 when 1 =>
      	    return X;
      	 when others =>
	    case Level is
	       when Value =>
		  return (Val  => X.Val ** K,
			  Grad => G0,
			  Hessian => H0);
	       when Gradient =>
		  Y := Real (K) * X.Val ** (K - 1);
		  return (Val  => X.Val ** K,
			  Grad => Y * X.Grad,
			  Hessian => H0);
	       when Hessian =>
		  Y := Real (K) * X.Val ** (K - 1);
		  Z := Real (K * (K - 1)) * X.Val ** (K - 2);
		  H := Y * X.Hessian + Z * Outer (X.Grad, X.Grad);
		  return (Val  => X.Val ** K,
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
	    return (Val => S,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (Val => S,
		    Grad => C * X.Grad,
		    Hessian => H0);
	 when Hessian =>
	    return (Val => S,
		    Grad => C * X.Grad,
		    Hessian => -S * Outer (X.Grad, X.Grad) + C * X.Hessian);
      end case;
   end Sin;
   function Cos (X : in AD_Type) return AD_Type is
      use Real_Functions;
      S : constant Real := Sin (X.Val);
      C : constant Real := Cos (X.Val);
   begin
      case Level is
	 when Value =>
	    return (Val  => C,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    return (Val  =>  C,
		    Grad => -S * X.Grad,
		    Hessian => H0);
	 when Hessian =>
	    return (Val  =>  C,
		    Grad => -S * X.Grad,
		    Hessian => -C * X.Grad * X.Grad - S * X.Hessian);
      end case;
   end Cos;
   function Tan (X : in AD_Type) return AD_Type is
      use Real_Functions;
      C : constant Real := Cos (X.Val);
      T : Real;
      Y : Real;
      Z : Real_Vector (1 .. N);
   begin
      pragma Assert (C /= 0.0);
      T := Sin (X.Val) / C;
      Y := (1.0 - T ** 2);
      case Level is
	 when Value =>
	    return (Val  => T,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    Z := Y * X.Grad;
	    return (Val  => T,
		    Grad => Z,
		    Hessian => H0);
	 when Hessian =>
	    Z := Y * X.Grad;
	    return (Val  => T,
		    Grad => Z,
		    Hessian => (1.0 + Y) * X.Hessian + ((2.0 * T) * Z) * X.Grad);
      end case;
   end Tan;
	      
      
   function Exp (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : constant Real := Exp (X.Val);
      G : Real_Vector (1 .. N);
   begin
      case Level is
	 when Value =>
	    return (Val  => Y,
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    G := Y * X.Grad;
	    return (Val  => Y,
		    Grad => G,
		    Hessian => H0);
	 when Hessian =>
	    G := Y * X.Grad;
	    return (Val  => Y,
		    Grad => G,
		    Hessian => Y * X.Hessian + G * X.Grad);
      end case;
   end Exp;
   
   
   function Log (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : Real;
      Z : Real_Vector (1 .. N);
   begin
      pragma Assert (X.Val > 0.0);
      Y := 1.0 / X.Val;
      case Level is
	 when Value =>
	    return (Val  => Log (X.Val),
		    Grad => G0,
		    Hessian => H0);
	 when Gradient =>
	    Z := Y * X.Grad;
	    return (Val  => Log (X.Val),
		    Grad => Z,
		    Hessian => H0);
	 when Hessian =>
	    Z := Y * X.Grad;
	    return (Val  => Log (X.Val),
		    Grad => Z,
		    Hessian => -Z * Z + Y * X.Hessian);
      end case; 
   end Log;
   
      
   function "-" (X : in AD_Type) return AD_Type is
   begin
      return (Val  => -X.Val,
   	      Grad => -X.Grad,
   	      Hessian => -X.Hessian);
   end "-";
   
   function "*" (X : in Real; Y : in AD_Type) return AD_Type is
   begin
      return (Val  => X * Y.Val,
   	      Grad => X * Y.Grad,
   	      Hessian => X * Y.Hessian);
   end "*";
   
   procedure print (X : in AD_Type) is
   begin
      null;
   end Print;
end Auto_Differentiation_Dense;
