all: develop latest v3.5.x v3.4.x v3.3.x v3.2.x

develop: develop-cpu develop-gpu

develop-cpu:
	docker build \
            -t systemsgenetics/kinc:develop-cpu \
            --build-arg NVIDIA_HEADLESS=1 \
            --build-arg ACE_REVISION=v3.2.0 \
            --build-arg KINC_REVISION=develop \
            --build-arg KINC_R_REVISION=v1.2 \
            .

develop-gpu:
	docker build \
            -t systemsgenetics/kinc:develop-gpu \
            --build-arg NVIDIA_HEADLESS=0 \
            --build-arg ACE_REVISION=v3.2.0 \
            --build-arg KINC_REVISION=develop \
            --build-arg KINC_R_REVISION=v1.2 \
            .


latest: latest-cpu latest-gpu

latest-cpu:
	docker build \
		-t systemsgenetics/kinc:latest-cpu \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=v3.2.0 \
		--build-arg KINC_REVISION=v3.4.2 \
		--build-arg KINC_R_REVISION=v1.2 \
		.

latest-gpu:
	docker build \
		-t systemsgenetics/kinc:latest-gpu \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=v3.2.0 \
		--build-arg KINC_REVISION=v3.4.2 \
		--build-arg KINC_R_REVISION=v1.2 \
		.
v3.5.x: v3.5.x-cpu v3.5.x-gpu

v3.5.x-cpu:
	docker build \
		-t systemsgenetics/kinc:3.5.0-cpu \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=v3.2.0 \
		--build-arg KINC_REVISION=v3.5.0 \
		--build-arg KINC_R_REVISION=v1.2 \
		.

v3.5.x-gpu:
	docker build \
		-t systemsgenetics/kinc:3.5.0-gpu \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=v3.2.0 \
		--build-arg KINC_REVISION=v3.5.0 \
		--build-arg KINC_R_REVISION=v1.2 \
		.


v3.4.x: v3.4.x-cpu v3.4.x-gpu

v3.4.x-cpu:
	docker build \
		-t systemsgenetics/kinc:3.4.2-cpu \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=v3.2.0 \
		--build-arg KINC_REVISION=v3.4.2 \
		--build-arg KINC_R_REVISION=v1.2 \
		.

v3.4.x-gpu:
	docker build \
		-t systemsgenetics/kinc:3.4.2-gpu \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=v3.2.0 \
		--build-arg KINC_REVISION=v3.4.2 \
		--build-arg KINC_R_REVISION=v1.2 \
		.

v3.3.x: v3.3.x-cpu v3.3.x-gpu

v3.3.x-cpu:
	docker build \
		-t systemsgenetics/kinc:3.3.0-cpu \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=v3.1.0 \
		--build-arg KINC_REVISION=v3.3.0 \
		.

v3.3.x-gpu:
	docker build \
		-t systemsgenetics/kinc:3.3.0-gpu \
		--build-arg NVIDIA_HEADLESS=0 \
		--build-arg ACE_REVISION=v3.1.0 \
		--build-arg KINC_REVISION=v3.3.0 \
		.

v3.2.x:
	docker build \
		-t systemsgenetics/kinc:3.2.2 \
		--build-arg NVIDIA_HEADLESS=1 \
		--build-arg ACE_REVISION=v3.0.2 \
		--build-arg KINC_REVISION=v3.2.2 \
		.

push:
	docker push systemsgenetics/kinc:develop-cpu
	docker push systemsgenetics/kinc:develop-gpu
	docker push systemsgenetics/kinc:latest-cpu
	docker push systemsgenetics/kinc:latest-gpu
	docker push systemsgenetics/kinc:3.4.2-cpu
	docker push systemsgenetics/kinc:3.4.2-gpu
	docker push systemsgenetics/kinc:3.3.0-cpu
	docker push systemsgenetics/kinc:3.3.0-gpu
	docker push systemsgenetics/kinc:3.2.2

pull:
	docker pull systemsgenetics/kinc:develop-cpu
	docker pull systemsgenetics/kinc:develop-gpu
	docker pull systemsgenetics/kinc:latest-cpu
	docker pull systemsgenetics/kinc:latest-gpu
	docker pull systemsgenetics/kinc:3.4.2-cpu
	docker pull systemsgenetics/kinc:3.4.2-gpu
	docker pull systemsgenetics/kinc:3.3.0-cpu
	docker pull systemsgenetics/kinc:3.3.0-gpu
	docker pull systemsgenetics/kinc:3.2.2

clean:
	docker image rm -f systemsgenetics/kinc:develop-cpu
	docker image rm -f systemsgenetics/kinc:develop-gpu
	docker image rm -f systemsgenetics/kinc:latest-cpu
	docker image rm -f systemsgenetics/kinc:latest-gpu
	docker image rm -f systemsgenetics/kinc:3.4.2-cpu
	docker image rm -f systemsgenetics/kinc:3.4.2-gpu
	docker image rm -f systemsgenetics/kinc:3.3.0-cpu
	docker image rm -f systemsgenetics/kinc:3.3.0-gpu
	docker image rm -f systemsgenetics/kinc:3.2.2
