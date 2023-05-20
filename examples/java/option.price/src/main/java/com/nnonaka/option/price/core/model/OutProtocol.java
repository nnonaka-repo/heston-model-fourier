package com.nnonaka.option.price.core.model;

import java.io.Serializable;
import java.util.List;

public class OutProtocol implements Serializable {
    private List<KeyValueData> header;
    private Object data;
    private int statusCode;

    public List<KeyValueData> getHeader() {
        return header;
    }

    public void setHeader(List<KeyValueData> header) {
        this.header = header;
    }

    public Object getData() {
        return data;
    }

    public void setData(Object data) {
        this.data = data;
    }

    public int getStatusCode() {
        return statusCode;
    }

    public void setStatusCode(int statusCode) {
        this.statusCode = statusCode;
    }
}
