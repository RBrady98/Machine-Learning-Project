APP_NAME="java_app"

build:
	docker build -t $(APP_NAME) .

run:
	docker run -it --rm $(APP_NAME)
