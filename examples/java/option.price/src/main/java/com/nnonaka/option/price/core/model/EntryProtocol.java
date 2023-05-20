package com.nnonaka.option.price.core.model;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.nnonaka.option.price.core.model.enumeration.Method;

import java.io.Serializable;
import java.util.List;

public class EntryProtocol implements Serializable {
    @JsonProperty("query")
    private List<KeyValueData> query;
    @JsonProperty("header")
    private List<KeyValueData> header;
    @JsonProperty("path")
    private List<KeyValueData> path;
    @JsonProperty("body")
    private Object body;
    @JsonProperty("method")
    private Method method;

    public List<KeyValueData> getQuery() {
        return query;
    }

    public void setQuery(List<KeyValueData> query) {
        this.query = query;
    }

    public List<KeyValueData> getHeader() {
        return header;
    }

    public void setHeader(List<KeyValueData> header) {
        this.header = header;
    }

    public List<KeyValueData> getPath() {
        return path;
    }

    public void setPath(List<KeyValueData> path) {
        this.path = path;
    }


    public Object getBody() {
        return body;
    }

    public void setBody(Object body) {
        this.body = body;
    }

    public Method getMethod() {
        return method;
    }

    public void setMethod(Method method) {
        this.method = method;
    }
}
