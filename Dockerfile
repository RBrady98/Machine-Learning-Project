FROM openjdk:12-alpine
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp
RUN javac main.java
CMD ["java", "main"]
