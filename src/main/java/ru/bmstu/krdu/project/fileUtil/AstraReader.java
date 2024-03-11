package ru.bmstu.krdu.project.fileUtil;

import ru.bmstu.krdu.project.dto.AstraResult;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class AstraReader {

    private static final ClassLoader classLoader = AstraReader.class.getClassLoader();
    private final Map<String, Pattern> patternMap = new HashMap<>();

    public AstraReader() {
        patternMap.put("p", Pattern.compile("P=(\\d+(?:\\.\\d+)?) ")); //p должно быть первым!!!!
        patternMap.put("k", Pattern.compile("k=(\\d+(?:\\.\\d+)?) "));
        patternMap.put("t", Pattern.compile("T=(\\d+(?:\\.\\d+)?) "));
        patternMap.put("r", Pattern.compile("R\\.\\?=(\\d+(?:\\.\\d+)?)"));
        patternMap.put("i", Pattern.compile("I\\?=(\\d+(?:\\.\\d+)?)"));
        patternMap.put("foth", Pattern.compile("F/F\\*=(\\d+(?:\\.\\d+)?)"));
        patternMap.put("fflow", Pattern.compile("F\"=(\\d+(?:\\.\\d+)?)"));
    }

    public List<AstraResult> readParams(String path) {
        LinkedList<AstraResult> result = new LinkedList<>();

        try(InputStream inputStream = classLoader.getResourceAsStream(path)) {
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream, StandardCharsets.ISO_8859_1));
            String line;
            while ((line = reader.readLine()) != null) {
                Map<String, Matcher> matchers = new LinkedHashMap<>();
                String finalLine = line;
                patternMap.forEach((key, value) -> matchers.put(key, value.matcher(finalLine)));
                for (Map.Entry<String, Matcher> es: matchers.entrySet()
                     ) {
                    AstraResult astraResult;
                    String param = es.getKey();
                    Matcher matcher = es.getValue();
                    if(param.equals("p") && matcher.find()) {
                        astraResult = new AstraResult();
                        astraResult.setField("p", Double.parseDouble(matcher.group(1)));
                        result.add(astraResult);
                    } else if(!result.isEmpty()){
                        astraResult = result.get(result.size()-1);
                    } else {
                        continue;
                    }

                    if(matcher.find()) {
                        astraResult.setField(param, Double.parseDouble(matcher.group(1)));
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return result;
    }
}
