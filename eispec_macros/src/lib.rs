extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, Data, DeriveInput, Fields, GenericArgument, PathArguments, Type};

///
/// Macro that automatically defines the data getter/setter methods required for `ImpedanceData<T>`.
///
#[proc_macro_attribute]
pub fn impl_impedance_data_accessors(_attr: TokenStream, item: TokenStream) -> TokenStream {
    let input = parse_macro_input!(item as DeriveInput);
    let struct_name = &input.ident;

    let expanded = match input.data {
        Data::Struct(ref data_struct) => match &data_struct.fields {
            Fields::Unnamed(fields) => {
                if fields.unnamed.len() == 1 {
                    let field_type = &fields.unnamed[0].ty;
                    let numeric_type = extract_inner_type(field_type);

                    if let Some(numeric_type) = numeric_type {
                        quote! {
                            #input

                            impl ImpedanceDataAccessors<#numeric_type> for #struct_name {
                                fn get_freqs(&self) -> &[Frequency<#numeric_type>] {
                                    self.0.get_freqs()
                                }

                                fn get_zmeas(&self) -> &[Impedance<#numeric_type>] {
                                    self.0.get_zmeas()
                                }

                                fn get_zerr(&self) -> &[Impedance<#numeric_type>] {
                                    self.0.get_zerr()
                                }

                                fn get_parameters(&self) -> Option<&[ModelParameter<#numeric_type>]> {
                                    self.0.get_parameters()
                                }

                                fn set_parameters(&mut self, parameters: Vec<ModelParameter<#numeric_type>>) {
                                    self.0.set_parameters(parameters)
                                }

                                fn from_csv(filename: &str) -> Result<#struct_name, ImpedanceDataError> {
                                    ImpedanceData::<#numeric_type>::from_csv(filename).map(#struct_name)
                                }
                            }
                        }
                    } else {
                        quote! {
                            compile_error!("Struct should have exactly one unnamed field of type ImpedanceData<T>.");
                        }
                    }
                } else {
                    quote! {
                        compile_error!("Struct should have exactly one unnamed field of type ImpedanceData<T>.");
                    }
                }
            }
            _ => quote! {
                compile_error!("Struct should have exactly one unnamed field of type ImpedanceData<T>.");
            },
        },
        _ => quote! {
            compile_error!("This macro only supports structs.");
        },
    };

    TokenStream::from(expanded)
}

///
/// Extracts the inner type `T` from `ImpedanceData<T>`.
///
fn extract_inner_type(ty: &Type) -> Option<Type> {
    let segment = match ty {
        Type::Path(t) => t.path.segments.last()?,
        _ => return None,
    };
    if segment.ident != "ImpedanceData" {
        return None;
    }
    let angle_bracketed_args = match &segment.arguments {
        PathArguments::AngleBracketed(args) => args,
        _ => return None,
    };

    match angle_bracketed_args.args.first()? {
        GenericArgument::Type(inner_type) => Some(inner_type.clone()),
        _ => None,
    }
}
