#' Achieve numerical conversion between units commonly used in pharmacokinetics.
#' such as: mg to ng, nmol to mol, ml to L, noml to mg, nmol/L to ng/mL.
#'
#' @title pk.unit.trans
#' @param value Vector or 1 column form data frame, Must be numeric.
#' @param from Vector or 1 column form data frame, Must be character.
#' @param to Vector or 1 column form data frame, Must be character.
#' @param mw molecular weight,Vector or 1 column form data frame, Must be numeric.
#' @export
#' @return Numerical Vector.
#' @import stringi maditr
#' @examples
#' unit_trans(1,"nmol/L","ng/mL",mw=150000) #return 150
#' unit_trans(1,"nM","ng/mL",mw=150000) #return 150
#' unit_trans(1,"nmol","mg",mw=150000) #return 0.15
#' unit_trans(1,"mg","ug") #return 1000
#' unit_trans(1,"mL","uL") #return 1000

# library(stringi)
# library(maditr)

## 定义单位转换函数
unit_trans <- function(value,from,to,mw=NULL) {

	# 单位转换
	Value_Form <- value
	Input_Raw <- from
	Output_Raw <- to
	RefUnitTable<-RefUnitTable
	mw <- mw

	# 1清理格式
  	## 1.1清理空格
	Input_Value <- gsub(' ', '', Input_Raw)
	Output_Value <- gsub(' ', '', Output_Raw)

  	##1.2 将非标单位字符替换为标准的单位字符
	#   library(stringi)
	Input_Value_new   <- stri_replace_all_regex(str= Input_Value   ,pattern= Non_standard_unit_transform$raw_unit, replacement= Non_standard_unit_transform$source_unit,vectorize_all = FALSE)
	Output_Value_new  <- stri_replace_all_regex(str= Output_Value  ,pattern= Non_standard_unit_transform$raw_unit, replacement= Non_standard_unit_transform$source_unit,vectorize_all = FALSE)
	if(Input_Value_new!=Input_Value){
		#单位'%s'被标准化为单位'%s'。
		message(sprintf("The unit '%s' has been normalized to the unit '%s'.(\u5355\u4f4d'%s'\u88ab\u6807\u51c6\u5316\u4e3a\u5355\u4f4d'%s'\u3002)",Input_Value,Input_Value_new,Input_Value,Input_Value_new))
		Input_Value <- Input_Value_new
	}
	if(Output_Value_new!=Output_Value){
		#单位'%s'被标准化为单位'%s'。
		message(sprintf("The unit '%s' has been normalized to the unit '%s'.(\u5355\u4f4d'%s'\u88ab\u6807\u51c6\u5316\u4e3a\u5355\u4f4d'%s'\u3002)",Output_Value,Output_Value_new,Output_Value,Output_Value_new))
		Output_Value <- Output_Value_new
	}

	# 2 使用运算符号拆分字符串，仅使用乘法符号和除法符号拆分
	Input_VSplit <- unlist(strsplit(
		x=unlist(strsplit(x=Input_Value,"/",fixed=T)),
	"*",fixed=T))
	Output_VSplit <- unlist(strsplit(
		x=unlist(strsplit(x=Output_Value,"/",fixed=T)),
	"*",fixed=T))

	#3 忽略掉括号字符"("、")"，提取基本单位
	Input_Units <- gsub(
		x=gsub(x=Input_VSplit ,pattern="(", replacement="",fixed=T)
	,pattern=")", replacement="",fixed=T)
	Output_Units <- gsub(
		x=gsub(x=Output_VSplit ,pattern="(", replacement="",fixed=T)
	,pattern=")", replacement="",fixed=T)

	## 4.1 检查拆分后的基本单位是否在单位参考表中有包含，如果未包含则提示当前不支持对该单位的转换
	if ( !all(Input_Units %in% RefUnitTable$source_unit)) {
		stop_t1 <- paste0(Input_Units[!(Input_Units %in% RefUnitTable$source_unit)],collapse = ",")
		#不支持单位%s
		message(sprintf("%s units is not support.(\u4e0d\u652f\u6301\u5355\u4f4d %s \u3002)",stop_t1,stop_t1))
		return() ##stop()函数会报奇怪的错误，所以这里改用message()，并返回空值
	}
	if ( !all(Output_Units %in% RefUnitTable$source_unit)) {
		stop_t1 <- paste0(Output_Units[!(Output_Units %in% RefUnitTable$source_unit)],collapse = ",")
		#不支持单位%s
		message(sprintf("%s units is not support.(\u4e0d\u652f\u6301\u5355\u4f4d %s \u3002)",stop_t1,stop_t1))
		return() ##stop()函数会报奇怪的错误，所以这里改用message()，并返回空值
	}

	## 4.2 检查拆分后的元素数量是否一致如果不一致则报错，提示无法转换
	F_length <- length(Input_Units) == length(Output_Units)
	if (!F_length) {
		#无法由%s单位转换为%s单位。
		message(sprintf("Cannot be converted from %s units to %s units.(\u65e0\u6cd5\u7531%s\u5355\u4f4d\u8f6c\u6362\u4e3a%s\u5355\u4f4d\u3002)",Input_Value,Output_Value,Input_Value,Output_Value))
		return() ##stop()函数会报奇怪的错误，所以这里改用message()，并返回空值
	}

	# 4.3 拼接 Input_Units 和 Output_Units 为 IO_Units
	IO_Units <-paste(Input_Units, Output_Units, sep ="_")

	#5 使用 IO_Units 从数据库中查查找换算系数，并将储存换算系数的向量命名为 Unit_ConversionFactors
	# library(maditr)
	Unit_ConversionFactors <-  vlookup(IO_Units, RefUnitTable ,19)
	Unit_bracket_ConversionFactors <-paste0("(",Unit_ConversionFactors, ")",sep ="")

	#6 在 Input_Value 中查找 Input_Units 并替换为 Unit_ConversionFactors ，然后将查找替换后的值储存在命名为 InputTransExprs 的变量中
	#library(stringi)
	InputTransExprs <- stri_replace_all_regex(str= Input_Value, pattern= Input_Units, replacement= Unit_bracket_ConversionFactors,vectorize_all = FALSE)

	## 7.1 检查转换表达式中是否需要分子量mw，如果需要则进一步检查是否输入了mw，以及mw是否是数值。
	F_mw <- grepl("mw", InputTransExprs)
	if (F_mw) {
		if (is.null(mw)){
			#(分子量(mw)当前为null，请输入一个数值。)
			message("Molecular Weight(mw) currently is NULL, please enter a number.(\u5206\u5b50\u91cf(mw)\u5f53\u524d\u4e3anull\uff0c\u8bf7\u8f93\u5165\u4e00\u4e2a\u6570\u503c\u3002)") 
			return()
		}
		if(!is.numeric(mw)) {
			#(分子量(mw)当前不是数值，请输入一个数值。)
			message("Molecular Weight(mw) currently not is number, please enter a number.(\u5206\u5b50\u91cf(mw)\u5f53\u524d\u4e0d\u662f\u6570\u503c\uff0c\u8bf7\u8f93\u5165\u4e00\u4e2a\u6570\u503c\u3002)")
			return()
		}
		if(mw<=0) {
			#分子量(mw)当前不是正数(<=0)，请输入一个正数。)
			message("Molecular Weight(mw) currently not is positive number(<=0), please enter a positive number.(\u5206\u5b50\u91cf(mw)\u5f53\u524d\u4e0d\u662f\u6b63\u6570(<=0)\uff0c\u8bf7\u8f93\u5165\u4e00\u4e2a\u6b63\u6570\u3002)")
			return()
		}
	}

	#8 value + InputTransExprs合并，然后求出数值
	Exprs <- paste0(Value_Form , "*" , InputTransExprs)
	Value_To <- eval(parse(text=Exprs))

	return(Value_To)
}
