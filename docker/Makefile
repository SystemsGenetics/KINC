all:
	docker build \
		-t systemsgenetics/kinc:latest-cpu \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=develop \
		--build-arg KINC_REVISION=master \
		.
	docker build \
		-t systemsgenetics/kinc:latest-gpu \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=develop \
		--build-arg KINC_REVISION=master \
		.
	docker build \
		-t systemsgenetics/kinc:3.3.0-cpu \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=v3.1.0 \
		--build-arg KINC_REVISION=v3.3.0 \
		.
	docker build \
		-t systemsgenetics/kinc:3.3.0-gpu \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=v3.1.0 \
		--build-arg KINC_REVISION=v3.3.0 \
		.
	docker build \
		-t systemsgenetics/kinc:3.2.2 \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=v3.0.2 \
		--build-arg KINC_REVISION=v3.2.2 \
		.

push:
	docker push systemsgenetics/kinc:latest-cpu
	docker push systemsgenetics/kinc:latest-gpu
	docker push systemsgenetics/kinc:3.3.0-cpu
	docker push systemsgenetics/kinc:3.3.0-gpu
	docker push systemsgenetics/kinc:3.2.2

pull:
	docker pull systemsgenetics/kinc:latest-cpu
	docker pull systemsgenetics/kinc:latest-gpu
	docker pull systemsgenetics/kinc:3.3.0-cpu
	docker pull systemsgenetics/kinc:3.3.0-gpu
	docker pull systemsgenetics/kinc:3.2.2

clean:
	docker image rm -f systemsgenetics/kinc:latest-cpu
	docker image rm -f systemsgenetics/kinc:latest-gpu
	docker image rm -f systemsgenetics/kinc:3.3.0-cpu
	docker image rm -f systemsgenetics/kinc:3.3.0-gpu
	docker image rm -f systemsgenetics/kinc:3.2.2
