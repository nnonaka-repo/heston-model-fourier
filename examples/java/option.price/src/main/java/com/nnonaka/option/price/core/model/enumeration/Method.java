package com.nnonaka.option.price.core.model.enumeration;

public enum Method {
    GET("GET"),
    POST("POST"),
    PUT("PUT"),
    DELETE("DELETE"),
    PATCH("PATCH"),
    HEAD("HEAD"),
    OPTIONS("OPTIONS"),
    TRACE("TRACE");

    private String tipo;

    private Method(String tipo){ this.tipo=tipo;}
    public static Method getMethod(String tipo){
        for(Method method: Method.values()){
            if(method.tipo.equals(tipo)){
                return method;
            }
        }
        return null;
    }
}
