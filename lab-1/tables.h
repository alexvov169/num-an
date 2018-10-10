#pragma once

#include <iostream>

template<typename ApproximatorType, typename ValueType, typename IntegerType>
class tables {
public:
	typedef ApproximatorType approximator_type;
	typedef ValueType value_type;
	typedef IntegerType integer_type;
	tables(approximator_type approximator,
		value_type second_table_error,
		value_type initial_error,
		value_type final_error,
		value_type error_step,
		value_type first_table_value,
		value_type initial_value,
		value_type final_value,
		integer_type nvalues) :
		second_table_error(second_table_error),
		initial_error(initial_error),
		final_error(final_error),
		error_step(error_step),
		first_table_value(first_table_value),
		initial_value(initial_value),
		final_value(final_value),
		nvalues(nvalues),
		approximator(approximator) {}
private:
	void print_first_table_row(std::ostream& stream, value_type error) const {
		value_type approximated = approximator(first_table_value, error); // must be called first
		stream << std::scientific << std::right
			<< std::setw(20) << error
			<< std::setw(20) << approximator.get_nterms()
			<< std::setw(20) << approximator.get_remainder_term()
			<< std::setw(20) << approximator.get_exact(first_table_value) - approximated
			<< std::endl;
	}
	void print_first_table(std::ostream& stream) const {
		stream
			<< std::setw(20) << "eps"
			<< std::setw(20) << "n"
			<< std::setw(20) << "remainder term"
			<< std::setw(20) << "delta"
			<< std::endl;
		auto mult_e_min3 = [this](const value_type& x) { return error_step * x; };
		for (const auto& error : rrange(initial_error, second_table_error,
			mult_e_min3, std::greater<double>())) {
			print_first_table_row(stream, error);
		}
		print_first_table_row(stream, second_table_error);
		second_table_nterms = approximator.get_nterms();
		for (const auto& error : rrange(second_table_error * error_step, final_error,
			mult_e_min3, std::greater<double>())) {
			print_first_table_row(stream, error);
		}
	}
	void print_second_table(std::ostream& stream) const {
		stream
			<< std::setw(20) << "x"
			<< std::setw(20) << "remainder term"
			<< std::setw(20) << "delta"
			<< std::endl;
		value_type approximated;
		for (const auto& x : range(initial_value, final_value, nvalues)) {
			approximated = approximator(first_table_value, second_table_nterms);
			stream << std::scientific << std::right
				<< std::setw(20) << x
				<< std::setw(20) << approximator.get_remainder_term()
				<< std::setw(20) << approximator.get_exact(x) - approximated
				<< std::setw(20) << approximator.get_exact(x)
				<< std::endl;
		}
	}
private:
	value_type second_table_error;
	value_type initial_error;
	value_type final_error;
	value_type error_step;
	value_type first_table_value;
	value_type initial_value;
	value_type final_value;
	integer_type nvalues;
	approximator_type approximator;

	mutable integer_type second_table_nterms;
	template<typename ApprxTy, typename ValTy, typename IntTy>
	friend std::ostream& operator<<(std::ostream&, const tables<ApprxTy, ValTy, IntTy>&);
};

template<typename ApprxTy, typename ValTy, typename IntTy>
inline auto make_tables(ApprxTy approximator,
	ValTy second_table_error,
	ValTy initial_error,
	ValTy final_error,
	ValTy error_step,
	ValTy first_table_value,
	ValTy initial_value,
	ValTy final_value,
	IntTy nvalues) {
	return tables<ApprxTy, ValTy, IntTy>(approximator,
		second_table_error, initial_error, final_error, error_step,
		first_table_value, initial_value, final_value, nvalues);
}

template<typename ApprxTy, typename ValTy, typename IntTy>
std::ostream& operator<<(std::ostream& stream, const tables<ApprxTy, ValTy, IntTy>& object) {
	object.print_first_table(stream);
	stream << std::endl;
	object.print_second_table(stream);
	return stream;
}
