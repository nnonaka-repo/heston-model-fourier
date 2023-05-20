package com.nnonaka.option.price.core.model;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.io.Serializable;

public class KeyValueData implements Serializable {
    @JsonProperty("key")
    public String key;
    @JsonProperty("value")
    public String value;

    public String getKey() {
        return key;
    }

    public void setKey(String key) {
        this.key = key;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }
}
